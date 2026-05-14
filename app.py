import streamlit as st
import json
import math

# --- 1. 데이터 로드 및 보간 함수 ---
def load_json(filename):
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            return json.load(f)
    except:
        return None

def interpolate(temp, properties, key):
    try:
        sorted_props = sorted(properties, key=lambda x: x['temp'])
        temps = [float(p['temp']) for p in sorted_props]
        if temp <= temps[0]: return sorted_props[0][key]
        if temp >= temps[-1]: return sorted_props[-1][key]
        for i in range(len(temps)-1):
            t1, t2 = temps[i], temps[i+1]
            if t1 <= temp <= t2:
                v1, v2 = sorted_props[i][key], sorted_props[i+1][key]
                return v1 + (v2 - v1) * (temp - t1) / (t2 - t1)
    except:
        return 0.0
    return 0.0
    
# --- 2. 경제적 최적 지름 계산 엔진 (시행착오법) ---
def solve_economic_diameter(rho, mu, m_dot, c1, c2, t, n, a, b, f_multiplier, eta, epsilon):
    D_guess = 0.04 
    tolerance = 0.00001 
    max_iter = 50
    for i in range(max_iter):
        v = (4 * m_dot) / (rho * math.pi * D_guess**2)
        re = (4 * m_dot) / (math.pi * D_guess * mu) 
        if re > 2300:
            term = (epsilon / D_guess / 3.7)**1.11 + (6.9 / re)
            f = (-1.8 * math.log10(term))**-2 
        else:
            f = 64 / re if re > 0 else 0.01
        
        numerator = 40 * f * (m_dot**3) * (c2 / 1000) * t 
        denominator = n * (a + b) * (1 + f_multiplier) * c1 * eta * (math.pi**2) * (rho**2)
        D_new = (numerator / denominator)**(1 / (n + 5))
        if abs(D_new - D_guess) < tolerance:
            return D_new, f, re
        D_guess = D_new 
    return D_guess, f, re

# --- 3. 페이지 설정 및 데이터 로드 ---
st.set_page_config(page_title="공학용 유체 설계 시스템 v8.8", layout="wide")
st.title("🚀 공학용 유체 수송 설계 시스템 (Web v8.8)")
st.markdown("### 이론적 최적화 도출 및 상용 규격/비용 분석 통합 시스템")

f_db = load_json('fluids_db.json')
p_db = load_json('pipe_db.json')

if not f_db or not p_db:
    st.error("데이터 파일(JSON)을 찾을 수 없습니다.")
    st.stop()

# --- 4. 사이드바 설정 ---
with st.sidebar:
    st.header("[1] 유체 물성 및 경제성 설정")
    fluid_name = st.selectbox("대상 유체 선택", [f['name'] for f in f_db['fluids']])
    fluid_data = next(item for item in f_db['fluids'] if item['name'] == fluid_name)
    temp_range = [float(p['temp']) for p in fluid_data['properties']]
    min_t, max_t = min(temp_range), max(temp_range)
    
    fix_temp = st.checkbox("상온 고정 (20.0°C)", value=True)
    target_temp = 20.0 if fix_temp else st.number_input(f"운전 온도 ({min_t}~{max_t}°C)", min_value=min_t, max_value=max_t, value=20.0 if min_t <= 20 <= max_t else min_t)

    rho = interpolate(target_temp, fluid_data['properties'], 'rho')
    mu = interpolate(target_temp, fluid_data['properties'], 'mu')
    st.info(f"**밀도:** {rho:.2f} kg/m³  \n**점도:** {mu:.6f} Pa·s")

    st.divider()
    cost_grade_list = [g['grade'] for g in p_db['cost_grades']]
    sel_grade = st.selectbox("비용 등급 (2026 CPI 보정)", cost_grade_list)
    grade_data = next(g for g in p_db['cost_grades'] if g['grade'] == sel_grade)
    c1_value = st.number_input("설치비 상수 (C1)", value=grade_data['c1'])
    n_exponent = st.number_input("비용 지수 (n)", value=grade_data['n'])
    
    c2 = st.number_input("에너지 비용 ($/kWh)", value=0.04) 
    t_year = st.number_input("연간 가동 시간 (hr/yr)", value=6000) 
    eff_pump = st.slider("펌프 효율 (η)", 0.1, 1.0, 0.75) 
    ann_a = st.number_input("자본상환율 (a)", value=0.143, format="%.3f") 
    ann_b = st.number_input("유지보수율 (b)", value=0.01) 
    cost_f = st.number_input("부속품 배수 (F)", value=7.0) 

# --- 5. 메인 입력 영역 ---
col1, col2 = st.columns(2)
with col1:
    st.subheader("[2] 배관 규격 및 수치")
    nps_list = [p['nps'] for p in p_db['pipe_standards']]
    sel_nps = st.selectbox("NPS 선택", nps_list)
    pipe_info = next(p for p in p_db['pipe_standards'] if p['nps'] == sel_nps)
    sel_sch = st.selectbox("Schedule 선택", list(pipe_info['schedules'].keys()))
    d_val = st.number_input("관 안지름(ID)", value=pipe_info['schedules'][sel_sch]['id'])
    d_unit = st.selectbox("직경 단위", ["mm", "m", "inch"])
    l_val = st.number_input("배관 직선 거리", value=10.0)
    l_unit = st.selectbox("거리 단위", ["m", "km"])

with col2:
    st.subheader("[3] 유동 파라미터")
    v_val = st.number_input("유속/유량 입력", value=1.0)
    v_unit = st.selectbox("단위", ["m/s", "m³/s", "L/min", "L/s"])
    n_elbow = st.number_input("엘보 개수", min_value=0, value=0)
    n_valve = st.number_input("밸브 개수", min_value=0, value=0)

# --- 6. 계산 및 결과 출력 ---
if st.button("🚀 설계 시뮬레이션 실행", use_container_width=True):
    # 단위 환산
    D_current = d_val/1000 if d_unit=="mm" else (d_val*0.0254 if d_unit=="inch" else d_val)
    L = l_val*1000 if l_unit=="km" else l_val
    area = math.pi * (D_current**2) / 4
    if v_unit == "m/s": v = v_val
    elif v_unit == "m³/s": v = v_val / area
    elif v_unit == "L/min": v = (v_val/60000)/area
    else: v = (v_val/1000)/area
    m_dot = rho * v * area 

    # [단계 1] 이론적 D_opt 계산
    d_opt_m, f_opt, re_opt = solve_economic_diameter(rho, mu, m_dot, c1_value, c2, t_year, n_exponent, ann_a, ann_b, cost_f, eff_pump, 0.000046)

    # [단계 2] 상용 배관 규격 추천 (DB 매칭)
    recommended_pipe = None
    min_diff = float('inf')
    for p in p_db['pipe_standards']:
        if sel_sch in p['schedules']:
            db_id = p['schedules'][sel_sch]['id'] / 1000
            diff = abs(db_id - d_opt_m)
            if diff < min_diff:
                min_diff = diff
                recommended_pipe = {"nps": p['nps'], "id": db_id}

    # [단계 3] 실제 추천 배관 기반 경제성 산출
    D_real = recommended_pipe['id']
    v_real = (4 * m_dot) / (rho * math.pi * D_real**2)
    re_real = (rho * v_real * D_real) / mu
    f_real = (-1.8 * math.log10((0.000046/D_real/3.7)**1.11 + (6.9/re_real)))**-2 if re_real > 2300 else 64/re_real
    dp = f_real * (L / D_real) * (rho * v_real**2 / 2)
    power_kw = (dp * (m_dot / rho)) / (eff_pump * 1000)
    op_cost = power_kw * t_year * c2
    fixed_cost = c1_value * (D_real**n_exponent) * (ann_a + ann_b) * (1 + cost_f) * L
    tac = op_cost + fixed_cost

    st.divider()

    # --- 화면 출력 1: 이론적 도출 과정 (마음에 들어 하셨던 부분) ---
    st.subheader("💰 1. 이론적 최적 지름($D_{opt}$) 도출 근거")
    st.success(f"수치 해석 결과, 경제적 최적 지름은 **{d_opt_m*1000:.2f} mm** 입니다.")
    with st.expander("🔍 시행착오법 수렴 리포트 및 수식", expanded=True):
        st.latex(rf"D_{{opt}} = \left[ \frac{{40 \cdot {f_opt:.4f} \cdot {m_dot:.2f}^3 \cdot {c2/1000:.6f} \cdot {t_year}}}{{{n_exponent} \cdot ({ann_a:.3f} + {ann_b:.2f}) \cdot (1 + {cost_f:.1f}) \cdot {c1_value} \cdot {eff_pump} \cdot \pi^2 \cdot {rho:.0f}^2}} \right]^{{\frac{{1}}{{{n_exponent}+5}}}}")
        st.write(f"- 수렴 Reynolds No: {re_opt:.1f} | 수렴 마찰계수(f): {f_opt:.4f}")

    # --- 화면 출력 2: 상용 추천 및 연간 비용 (교수님 피드백 반영) ---
    st.subheader("📋 2. 상용 규격 권고 및 연간 총 비용(TAC)")
    res1, res2 = st.columns(2)
    res1.info(f"**추천 규격:** NPS {recommended_pipe['nps']} (Sch.{sel_sch})  \n- 실지름: {D_real*1000:.2f} mm")
    res2.metric("총 연간 비용 (TAC)", f"$ {tac:,.2f} /yr")
    
    st.write("### 연간 소요 비용 분석")
    st.bar_chart({"에너지비(OPEX)": op_cost, "설비비(CAPEX 환산)": fixed_cost})

    # 유동 상태 표시
    if 2300 < re_real < 4000: st.warning(f"⚠️ 현재 추천 배관 운전 시 **천이 영역**에 해당합니다 (Re={re_real:.1f})")
    st.info(f"설계 유속: {v_real:.3f} m/s | 흐름 상태: {'난류' if re_real>4000 else '층류' if re_real<=2300 else '천이'}")
