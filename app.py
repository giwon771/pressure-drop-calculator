import streamlit as st
import json
import math

# --- 데이터 로드 함수 ---
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
                v1 = sorted_props[i][key]
                v2 = sorted_props[i+1][key]
                return v1 + (v2 - v1) * (temp - t1) / (t2 - t1)
    except:
        return 0.0
    return 0.0
    
# --- 경제적 최적 지름 계산 엔진 (시행착오법) ---
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
            f = 64 / re
        
        numerator = 40 * f * (m_dot**3) * (c2 / 1000) * t 
        denominator = n * (a + b) * (1 + f_multiplier) * c1 * eta * (math.pi**2) * (rho**2)
        D_new = (numerator / denominator)**(1 / (n + 5))
        
        if abs(D_new - D_guess) < tolerance:
            return D_new, f, re
        D_guess = D_new 
    return D_guess, f, re

# --- 페이지 설정 ---
st.set_page_config(page_title="공학용 유체 설계 시스템 v8.7", layout="wide")
st.title("🚀 공학용 유체 수송 설계 시스템 (Web v8.7)")
st.markdown("### 상용 배관 규격 추천 및 연간 총 비용(TAC) 분석 시스템")

f_db = load_json('fluids_db.json')
p_db = load_json('pipe_db.json')

if not f_db or not p_db:
    st.error("데이터 파일(JSON)을 찾을 수 없습니다.")
    st.stop()

# --- 사이드바 설정 ---
with st.sidebar:
    st.header("[1] 유체 물성 설정")
    fluid_name = st.selectbox("대상 유체 선택", [f['name'] for f in f_db['fluids']])
    
    fluid_data = next(item for item in f_db['fluids'] if item['name'] == fluid_name)
    temp_range = [float(p['temp']) for p in fluid_data['properties']]
    min_t, max_t = min(temp_range), max(temp_range)
    
    fix_temp = st.checkbox(f"상온 고정 (20.0°C)", value=True)
    
    if fix_temp:
        target_temp = 20.0
        target_temp = max(min_t, min(target_temp, max_t))
    else:
        target_temp = st.number_input(
            f"운전 온도 ({min_t}°C ~ {max_t}°C)", 
            min_value=min_t, 
            max_value=max_t, 
            value=20.0 if min_t <= 20.0 <= max_t else min_t,
            help="물리적 데이터가 존재하여 신뢰성이 확보된 온도 범위입니다."
        )

    rho = interpolate(target_temp, fluid_data['properties'], 'rho')
    mu = interpolate(target_temp, fluid_data['properties'], 'mu')
    st.info(f"**밀도:** {rho:.2f} kg/m³\n\n**점도:** {mu:.6f} Pa·s")

    st.divider()
    st.header("[4] 경제성 분석 설정")
    
    cost_grade_list = [g['grade'] for g in p_db['cost_grades']]
    sel_grade = st.selectbox("배관 비용 등급 선택 (2026 CPI 보정)", cost_grade_list)
    grade_data = next(g for g in p_db['cost_grades'] if g['grade'] == sel_grade)
    
    c1_value = st.number_input("배관 비용 상수 (C1)", value=grade_data['c1'])
    n_exponent = st.number_input("비용 지수 (n)", value=grade_data['n'])
    
    st.markdown("---")
    c2 = st.number_input("에너지 비용 (C2, $/kWh)", value=0.04) 
    t_year = st.number_input("연간 가동 시간 (t, hr/yr)", value=6000) 
    eff_pump = st.slider("펌프 효율 (η)", 0.1, 1.0, 0.75) 
    
    col_eco1, col_eco2 = st.columns(2)
    with col_eco1:
        ann_a = st.number_input("자본상환율 (a)", value=0.143, format="%.3f") 
        ann_b = st.number_input("유지보수율 (b)", value=0.01) 
    with col_eco2:
        cost_f = st.number_input("부속품 배수 (F)", value=7.0) 

# --- 메인 입력창 ---
col1, col2 = st.columns(2)
with col1:
    st.subheader("[2] 배관 규격 및 수치")
    nps_list = [p['nps'] for p in p_db['pipe_standards']]
    sel_nps = st.selectbox("NPS 선택", nps_list)
    
    pipe_info = next(p for p in p_db['pipe_standards'] if p['nps'] == sel_nps)
    sch_list = list(pipe_info['schedules'].keys())
    sel_sch = st.selectbox("Schedule 선택", sch_list)
    
    d_val = st.number_input("관 안지름(ID)", value=pipe_info['schedules'][sel_sch]['id'])
    d_unit = st.selectbox("직경 단위", ["mm", "m", "inch"])
    l_val = st.number_input("배관 직선 거리", value=10.0)
    l_unit = st.selectbox("거리 단위", ["m", "km"])

with col2:
    st.subheader("[3] 유동 및 부속품")
    v_val = st.number_input("유속 또는 유량 입력", value=1.0)
    v_unit = st.selectbox("입력 단위 선택", ["m/s", "m³/s", "L/min", "L/s"])
    n_elbow = st.number_input("엘보(90°) 개수", min_value=0, value=0)
    n_valve = st.number_input("게이트밸브 개수", min_value=0, value=0)

# --- 계산 및 결과 ---
if st.button("🚀 설계 분석 및 계산 실행", use_container_width=True):
    D_input = d_val / 1000 if d_unit == "mm" else (d_val * 0.0254 if d_unit == "inch" else d_val)
    L = l_val * 1000 if l_unit == "km" else l_val
    area = math.pi * (D_input**2) / 4

    if v_unit == "m/s": v = v_val
    elif v_unit == "m³/s": v = v_val / area
    elif v_unit == "L/min": v = (v_val / 60000) / area
    else: v = (v_val / 1000) / area
    
    m_dot = rho * v * area 

    # 이론적 최적 지름 계산
    d_opt_m, final_f, final_re = solve_economic_diameter(
        rho, mu, m_dot, c1_value, c2, t_year, n_exponent, ann_a, ann_b, cost_f, eff_pump, 0.000046
    )

    # 1. 상용 배관 규격 추천 로직 (Database Matching)
    recommended_pipe = None
    min_diff = float('inf')
    
    for p in p_db['pipe_standards']:
        if sel_sch in p['schedules']:
            db_id = p['schedules'][sel_sch]['id'] / 1000 # m 환산
            diff = abs(db_id - d_opt_m)
            if diff < min_diff:
                min_diff = diff
                recommended_pipe = {"nps": p['nps'], "id": db_id}

    # 2. 추천 배관 기반 경제성 재계산
    D_real = recommended_pipe['id']
    v_real = (4 * m_dot) / (rho * math.pi * D_real**2)
    re_real = (rho * v_real * D_real) / mu
    
    # 실제 마찰계수 산출
    term = (0.000046 / D_real / 3.7)**1.11 + (6.9 / re_real)
    f_real = (-1.8 * math.log10(term))**-2 if re_real > 2300 else 64/re_real
    
    # 연간 운영 비용 (Operating Cost)
    delta_p = f_real * (L / D_real) * (rho * v_real**2 / 2)
    power_kw = (delta_p * (m_dot / rho)) / (eff_pump * 1000)
    annual_energy_cost = power_kw * t_year * c2
    
    # 연간 고정 비용 (Fixed Cost)
    annual_fixed_cost = c1_value * (D_real**n_exponent) * (ann_a + ann_b) * (1 + cost_f) * L
    total_annual_cost = annual_energy_cost + annual_fixed_cost

    st.divider()
    
    # 결과 출력 섹션
    st.subheader("💰 최종 설계 권고 및 경제성 평가")
    st.success(f"이론적 최적 지름은 **{d_opt_m*1000:.2f} mm** 이며, DB 기반 추천 규격은 **NPS {recommended_pipe['nps']} (Sch.{sel_sch})** 입니다.")
    
    c_res1, c_res2, c_res3 = st.columns(3)
    c_res1.metric("추천 배관 실지름", f"{D_real*1000:.2f} mm")
    c_res2.metric("총 연간 비용 (TAC)", f"$ {total_annual_cost:,.2f}")
    c_res3.metric("설계 유속 (v)", f"{v_real:.3f} m/s")

    # 비용 비중 시각화
    st.write("### 연간 비용 구성 분석")
    cost_data = {"에너지 비용 (OPEX)": annual_energy_cost, "설치 비용 (CAPEX 환산)": annual_fixed_cost}
    st.bar_chart(cost_data)

    with st.expander("🔍 경제성 분석 상세 리포트"):
        st.markdown(f"""
        ### 1. 상용 규격 추천 근거
        - 이론적 최적 지름($D_{{opt}}$)인 {d_opt_m*1000:.2f} mm에 가장 근접한 규격을 데이터베이스에서 매칭했습니다.
        - 매칭 규격: NPS {recommended_pipe['nps']} (안지름 {D_real*1000:.2f} mm)
        
        ### 2. 총 연간 비용 (Total Annual Cost) 산출식
        $$TAC = C_{{energy}} + C_{{fixed}}$$
        - **운영비 ($C_{{energy}}$):** 실제 배관에서의 압력 강하를 기반으로 한 펌프 동력 비용
        - **고정비 ($C_{{fixed}}$):** 2026 CPI가 반영된 설치 비용의 연간 자본 환산액
        """)
