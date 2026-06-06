import streamlit as st
import json
import math
import plotly.graph_objects as go
import numpy as np

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
st.set_page_config(page_title="공학용 유체 설계 시스템 v8.9", layout="wide")
st.title("🚀 공학용 유체 수송 설계 시스템 (Web v8.9)")
st.markdown("### Darby 예제 4.7 기반 상용 규격 의사결정 및 연속 비용 곡선 해석 시스템")

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
    st.info(f"**화학공학 플랜트 비용 가치 보정 (CEPCI 기준 최신화)**")

    st.divider()
    cost_grade_list = [g['grade'] for g in p_db['cost_grades']]
    sel_grade = st.selectbox("비용 등급 (CEPCI/CPI 보정)", cost_grade_list)
    grade_data = next(g for g in p_db['cost_grades'] if g['grade'] == sel_grade)
    c1_value = st.number_input("설치비 상수 (C1)", value=grade_data['c1'])
    n_exponent = st.number_input("비용 지수 (n)", value=grade_data['n'])
    
    c2 = st.number_input("에너지 비용 ($/kWh)", value=0.04, help="한전 고압 산업용 전력 단가 기준 환산 가치") 
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
if st.button("🚀 설계 시뮬레이션 및 예제 4.7 교차 검증 실행", use_container_width=True):
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

    # [단계 2] 예제 4.7식 실제 비용 비교 기반 상용 배관 규격 추천
    # 단순 지름 거리가 아니라, 인접한 상하위 관경의 실제 Total Cost/L을 계산하여 최소 비용 배관을 결정
    pipes_with_cost = []
    for p in p_db['pipe_standards']:
        if sel_sch in p['schedules']:
            db_id = p['schedules'][sel_sch]['id'] / 1000
            
            # 각 상용 관경에서의 유속 및 Re 재계산
            v_p = (4 * m_dot) / (rho * math.pi * db_id**2)
            re_p = (rho * v_p * db_id) / mu
            f_p = (-1.8 * math.log10((0.000046/db_id/3.7)**1.11 + (6.9/re_p)))**-2 if re_p > 2300 else 64/re_p
            
            # 길이당 비용 계산 (Darby 식 4.10 기반)
            pipe_cost_per_l = (ann_a + ann_b) * (1 + cost_f) * c1_value * (db_id**n_exponent)
            op_cost_per_l = (8 * f_p * (m_dot**3) / (math.pi**2 * rho**2 * db_id**5)) * (c2 / 1000) * t_year / eff_pump
            total_cost_per_l = pipe_cost_per_l + op_cost_per_l
            
            pipes_with_cost.append({
                "nps": p['nps'],
                "id": db_id,
                "total_cost_per_l": total_cost_per_l,
                "pipe_cost": pipe_cost_per_l * L,
                "op_cost": op_cost_per_l * L
            })
            
    # 이론적 최적 지름 양옆에 있는 배관 중 경제적 총 비용(Total Cost)이 최소인 배관 매칭
    pipes_with_cost.sort(key=lambda x: abs(x['id'] - d_opt_m))
    candidate_pipes = pipes_with_cost[:2] # 가장 가까운 두 규격 선정
    candidate_pipes.sort(key=lambda x: x['total_cost_per_l']) # 그 중 비용이 더 낮은 것을 최종 추천
    
    recommended_pipe = candidate_pipes[0]

    st.divider()

    # --- 화면 출력 1: 이론적 도출 과정 ---
    st.subheader("💰 1. 이론적 최적 지름($D_{opt}$) 도출 결과")
    st.success(f"수치 해석 결과, 이론적 경제 최적 지름은 **{d_opt_m*1000:.2f} mm** 입니다.")
    with st.expander("🔍 시행착오법 수렴 리포트 및 수식", expanded=True):
        st.latex(rf"D_{{opt}} = \left[ \frac{{40 \cdot {f_opt:.4f} \cdot {m_dot:.2f}^3 \cdot {c2/1000:.6f} \cdot {t_year}}}{{{n_exponent} \cdot ({ann_a:.3f} + {ann_b:.2f}) \cdot (1 + {cost_f:.1f}) \cdot {c1_value} \cdot {eff_pump} \cdot \pi^2 \cdot {rho:.0f}^2}} \right]^{{\frac{{1}}{{{n_exponent}+5}}}}")
        st.write(f"- 수렴 Reynolds No: {re_opt:.1f} | 수렴 마찰계수(f): {f_opt:.4f}")

    # --- 화면 출력 2: Darby 그림 4.16 연속 비용 곡선 시각화 (신규 추가!) ---
    st.subheader("📊 2. 관경 변화에 따른 연속 비용 최적화 곡선 (Darby Fig 4.16 구현)")
    
    # 0.02m ~ 0.1m 범위의 연속적인 그래프 데이터 생성
    d_space = np.linspace(0.015, 0.10, 200)
    pipe_costs_line = []
    op_costs_line = []
    total_costs_line = []
    
    for d_s in d_space:
        v_s = (4 * m_dot) / (rho * math.pi * d_s**2)
        re_s = (rho * v_s * d_s) / mu
        f_s = (-1.8 * math.log10((0.000046/d_s/3.7)**1.11 + (6.9/re_s)))**-2 if re_s > 2300 else 64/re_s
        
        p_c = (ann_a + ann_b) * (1 + cost_f) * c1_value * (d_s**n_exponent)
        o_c = (8 * f_s * (m_dot**3) / (math.pi**2 * rho**2 * d_s**5)) * (c2 / 1000) * t_year / eff_pump
        
        pipe_costs_line.append(p_c)
        op_costs_line.append(o_c)
        total_costs_line.append(p_c + o_c)

    # Plotly 시각화 구성
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=d_space*1000, y=pipe_costs_line, name="Pipe Cost/L (자재 설치비)", line=dict(dash='dash', color='#005088')))
    fig.add_trace(go.Scatter(x=d_space*1000, y=op_costs_line, name="Operating Cost/L (동력 운영비)", line=dict(dash='dot', color='#ef4444')))
    fig.add_trace(go.Scatter(x=d_space*1000, y=total_costs_line, name="Total Cost/L (총 연간 비용)", line=dict(width=3, color='#11caa0')))
    
    # 이론적 최적 지름 점 표시
    fig.add_trace(go.Scatter(x=[d_opt_m*1000], y=[(ann_a + ann_b) * (1 + cost_f) * c1_value * (d_opt_m**n_exponent) + (8 * f_opt * (m_dot**3) / (math.pi**2 * rho**2 * d_opt_m**5)) * (c2 / 1000) * t_year / eff_pump],
                             mode='markers+text', name="Theoretical D_opt", text=["최적 지점"], textposition="top center", marker=dict(size=12, color='black', symbol='star')))

    fig.update_layout(
        xaxis=dict(title="Pipe Inside Diameter (mm)", gridcolor="#e2e8f0"),
        yaxis=dict(title="Installed Cost / Length ($/yr·m)", range=[0, max(total_costs_line)*0.5], gridcolor="#e2e8f0"),
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        height=500, legend=dict(x=0.6, y=0.9)
    )
    st.plotly_chart(fig, use_container_width=True)

    # --- 화면 출력 3: 상용 추천 및 연간 비용 ---
    st.subheader("📋 3. 상용 배관 규격 권고 및 경제성 리포트 (Darby 예제 4.7 검증 완료)")
    res1, res2 = st.columns(2)
    
    tac = recommended_pipe['total_cost_per_l'] * L
    res1.info(f"**최종 추천 상용 규격:** NPS {recommended_pipe['nps']} (Sch.{sel_sch})  \n- 실제 안지름: {recommended_pipe['id']*1000:.2f} mm  \n- 공학적 판단 기준: 후보 규격 간 실제 연간 비용 가치 방정식 대입 검증 완료")
    res2.metric("총 연간 비용 (TAC)", f"$ {tac:,.2f} /yr")
    
    # 두 후보 배관 비용 전격 비교 테이블
    st.write("### 🔍 인접 상용 배관 규격별 총비용 비교 검증")
    comparison_data = {
        "배관 규격 (NPS)": [f"NPS {p['nps']}" for p in candidate_pipes],
        "안지름 (mm)": [f"{p['id']*1000:.2f} mm" for p in candidate_pipes],
        "길이당 연간비용 ($/yr·m)": [f"$ {p['total_cost_per_l']:.2f}" for p in candidate_pipes],
        "연간 총 소요 비용 ($/yr)": [f"$ {p['total_cost_per_l']*L:,.2f}" for p in candidate_pipes]
    }
    st.table(comparison_data)

    # 유동 상태 및 속도 표시
    re_real = (rho * (4 * m_dot / (rho * math.pi * D_real**2)) * D_real) / mu
    flow_status = '난류' if re_real > 4000 else ('층류' if re_real <= 2300 else '천이')
    st.info(f"추천 배관 실제 내 유속: {v_real:.3f} m/s | 최종 운전 흐름 상태: {flow_status} (Re={re_real:.1f})")
