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
    
# --- 경제적 최적 지름 계산 엔진 (예제 4.6 로직) ---
def solve_economic_diameter(rho, mu, m_dot, c1, c2, t, n, a, b, f_multiplier, eta, epsilon):
    D_guess = 0.04 # 초기 가정값: 4cm [cite: 298]
    tolerance = 0.00001 
    max_iter = 50
    for i in range(max_iter):
        v = (4 * m_dot) / (rho * math.pi * D_guess**2)
        re = (4 * m_dot) / (math.pi * D_guess * mu) # Re 정의 [cite: 293, 294]
        if re > 2300:
            term = (epsilon / D_guess / 3.7)**1.11 + (6.9 / re)
            f = (-1.8 * math.log10(term))**-2 # 마찰계수 f 산출 [cite: 302]
        else:
            f = 64 / re
        # D_opt 결정 방정식 [cite: 212, 285, 287]
        numerator = 40 * f * (m_dot**3) * (c2 / 1000) * t 
        denominator = n * (a + b) * (1 + f_multiplier) * c1 * eta * (math.pi**2) * (rho**2)
        D_new = (numerator / denominator)**(1 / (n + 5))
        if abs(D_new - D_guess) < tolerance:
            return D_new, f, re
        D_guess = D_new 
    return D_guess, f, re

# --- 페이지 설정 ---
st.set_page_config(page_title="공학용 유체 설계 시스템 v8.0", layout="wide")
st.title("🚀 공학용 유체 수송 설계 시스템 (Web v8.0)")
st.markdown("### 배관 내 유동 분석 및 압력 손실 계산기")

f_db = load_json('fluids_db.json')
p_db = load_json('pipe_db.json')

if not f_db or not p_db:
    st.error("데이터 파일(JSON)을 찾을 수 없습니다.")
    st.stop()

# --- 사이드바 설정 ---
with st.sidebar:
    st.header("[1] 유체 물성 설정")
    fluid_name = st.selectbox("대상 유체 선택", [f['name'] for f in f_db['fluids']])
    fix_temp = st.checkbox("상온 고정 (20°C)", value=True)
    target_temp = 20.0 if fix_temp else st.number_input("운전 온도 (°C)", value=20.0)

    fluid_data = next(item for item in f_db['fluids'] if item['name'] == fluid_name)
    rho = interpolate(target_temp, fluid_data['properties'], 'rho')
    mu = interpolate(target_temp, fluid_data['properties'], 'mu')
    st.info(f"밀도: {rho:.2f} kg/m³\n\n점도: {mu:.6f} Pa·s")

    # [4] 경제성 분석 설정 섹션 추가 (이 부분이 화면에 나와야 합니다!)
    st.divider()
    st.header("[4] 경제성 분석 설정")
    c2 = st.number_input("에너지 비용 (C2, $/kWh)", value=0.04) # 예제 4.6 [cite: 262]
    t_year = st.number_input("연간 가동 시간 (t, hr/yr)", value=6000) # [cite: 264]
    eff_pump = st.slider("펌프 효율 (η)", 0.1, 1.0, 0.75) # [cite: 268]
    
    col_eco1, col_eco2 = st.columns(2)
    with col_eco1:
        ann_a = st.number_input("자본상환율 (a)", value=0.143, format="%.3f") # 1/7 [cite: 267]
        ann_b = st.number_input("유지보수율 (b)", value=0.01) # [cite: 267]
    with col_eco2:
        cost_f = st.number_input("부속품 배수 (F)", value=7.0) # [cite: 265]
        n_exponent = st.selectbox("비용 지수 (n)", [1.0, 1.2, 1.4], index=1) # [cite: 266]
    c1_value = st.number_input("배관 비용 상수 (C1)", value=700.0) # [cite: 263]

# --- 메인 입력창 ---
col1, col2 = st.columns(2)
with col1:
    st.subheader("[2] 배관 규격 및 수치")
    nps_list = [p['nps'] for p in p_db['pipe_standards']]
    sel_nps = st.selectbox("NPS 선택", nps_list)
    d_val = st.number_input("관 안지름(ID)", value=52.5)
    d_unit = st.selectbox("직경 단위", ["mm", "m", "inch"])
    l_val = st.number_input("배관 직선 거리(L)", value=10.0)
    l_unit = st.selectbox("거리 단위", ["m", "km"])

with col2:
    st.subheader("[3] 유동 및 부속품")
    v_val = st.number_input("유속 또는 유량 입력", value=1.0)
    v_unit = st.selectbox("입력 단위 선택", ["m/s", "m³/s", "L/min", "L/s"])
    n_elbow = st.number_input("엘보(90°) 개수", min_value=0, value=0)
    n_valve = st.number_input("게이트밸브 개수", min_value=0, value=0)

# --- 계산 및 결과 ---
if st.button("🚀 설계 분석 및 계산 실행", use_container_width=True):
    # 단위 환산
    D = d_val / 1000 if d_unit == "mm" else (d_val * 0.0254 if d_unit == "inch" else d_val)
    L = l_val * 1000 if l_unit == "km" else l_val
    area = math.pi * (D**2) / 4

    if v_unit == "m/s": v = v_val
    elif v_unit == "m³/s": v = v_val / area
    elif v_unit == "L/min": v = (v_val / 60000) / area
    else: v = (v_val / 1000) / area
    
    m_dot = rho * v * area # 질량유량 계산 

    # Re 및 흐름 판정
    re = (rho * v * D) / mu
    if re <= 2300: flow_type, f, color = "층류", 64/re, "blue"
    else: flow_type, f, color = "난류", 0.316*(re**-0.25), "red"

    # 경제적 최적 지름 계산 호출 [cite: 212]
    d_opt_m, final_f, final_re = solve_economic_diameter(
        rho, mu, m_dot, c1_value, c2, t_year, n_exponent, ann_a, ann_b, cost_f, eff_pump, 0.000046
    )
    
    # 결과 출력
    st.divider()
    st.subheader("💰 경제적 최적 지름(D_opt) 분석")
    st.success(f"최적 경제적 지름은 **{d_opt_m*1000:.2f} mm** 입니다.")
    st.info("선정 근거: 초기 설치비와 연간 에너지 비용의 합계 최소화 [cite: 14, 203]")
    
    c1, c2 = st.columns(2)
    c1.metric("현재 설정 Re", f"{re:.1f}")
    c2.metric("최적 설계 Re", f"{final_re:.1f}")
