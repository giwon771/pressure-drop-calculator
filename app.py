import streamlit as st
import json
import math
import pandas as pd

# --- 데이터 로드 함수 ---
def load_json(filename):
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            return json.load(f)
    except:
        return None

def interpolate(temp, properties, key):
    try:
        # 1. 만약 properties가 리스트 안에 들어있다면 첫 번째 요소를 사용
        if isinstance(properties, list):
            properties = properties[0]
            
        # 2. 온도 키값들을 가져와서 숫자 리스트로 변환
        # (문자열 "20"이든 숫자 20이든 모두 float으로 변환)
        raw_keys = list(properties.keys())
        temps = sorted([float(k) for k in raw_keys])
        
        # 3. 입력 온도가 범위를 벗어날 경우 최댓값/최솟값 반환
        if temp <= temps[0]:
            target_key = next(k for k in raw_keys if float(k) == temps[0])
            return properties[target_key][key]
        if temp >= temps[-1]:
            target_key = next(k for k in raw_keys if float(k) == temps[-1])
            return properties[target_key][key]
        
        # 4. 선형 보간 수행
        for i in range(len(temps) - 1):
            t1, t2 = temps[i], temps[i+1]
            if t1 <= temp <= t2:
                # 실제 딕셔너리에서 해당 숫자와 일치하는 키를 다시 찾음
                k1 = next(k for k in raw_keys if float(k) == t1)
                k2 = next(k for k in raw_keys if float(k) == t2)
                
                v1 = properties[k1][key]
                v2 = properties[k2][key]
                
                return v1 + (v2 - v1) * (temp - t1) / (t2 - t1)
                
    except Exception as e:
        # 에러 발생 시 디버깅을 위해 0을 반환하거나 로그 출력
        return 0.001 # 아주 작은 값 반환
    return 0

# --- 페이지 설정 ---
st.set_page_config(page_title="공학용 유체 설계 시스템 v8.0", layout="wide")

st.title("🚀 공학용 유체 수송 설계 시스템 (Web v8.0)")
st.markdown("### 배관 내 유동 분석 및 압력 손실 계산기")

# 데이터 로드
f_db = load_json('fluids_db.json')
p_db = load_json('pipe_db.json')

if not f_db or not p_db:
    st.error("데이터 파일(JSON)을 찾을 수 없습니다.")
    st.stop()

# --- 사이드바: 유체 및 온도 설정 [1] ---
with st.sidebar:
    st.header("[1] 유체 물성 설정")
    fluid_name = st.selectbox("대상 유체 선택", [f['name'] for f in f_db['fluids']])
    
    fix_temp = st.checkbox("상온 고정 (20°C)", value=True)
    if fix_temp:
        target_temp = 20.0
    else:
        target_temp = st.number_input("운전 온도 (°C)", value=20.0)

    fluid_data = next(item for item in f_db['fluids'] if item['name'] == fluid_name)
    rho = interpolate(target_temp, fluid_data['properties'], 'rho')
    mu = interpolate(target_temp, fluid_data['properties'], 'mu')
    
    st.info(f"선택 유체: {fluid_name}\n\n밀도: {rho:.2f} kg/m³\n\n점도: {mu:.6f} Pa·s")

# --- 메인 화면: 입력창 구성 [2] & [3] ---
col1, col2 = st.columns(2)

with col1:
    st.subheader("[2] 배관 규격 및 수치")
    
    # 배관 규격 선택
    nps_list = [p['nps'] for p in p_db['pipe_standards']]
    sel_nps = st.selectbox("NPS 선택", nps_list)
    sel_sch = st.selectbox("Schedule 선택", ["Sch 40", "Sch 80", "Sch 160"])
    
    # 자동 안지름 계산 (선택 시 기본값 제안)
    # 실제 구현 시에는 p_db에서 매칭되는 ID를 가져오는 로직 필요
    d_val = st.number_input("관 안지름(ID)", value=52.5, format="%.2f")
    d_unit = st.selectbox("직경 단위", ["mm", "m", "inch"], index=0)
    
    # 거리 입력
    l_val = st.number_input("배관 직선 거리(L)", value=10.0)
    l_unit = st.selectbox("거리 단위", ["m", "km"], index=0)

with col2:
    st.subheader("[3] 유동 및 부속품")
    
    # 유량/유속 입력
    v_val = st.number_input("유속 또는 유량 입력", value=1.0)
    v_unit = st.selectbox("입력 단위 선택", ["m/s", "m³/s", "L/min", "L/s", "BTU/hr (Water)"])
    
    st.markdown("---")
    # 부속품
    n_elbow = st.number_input("엘보(90°) 개수", min_value=0, value=0)
    n_valve = st.number_input("게이트밸브 개수", min_value=0, value=0)

# --- 계산 실행 버튼 ---
if st.button("🚀 설계 분석 및 계산 실행", use_container_width=True):
    # [1. 단위 변환]
    D = d_val / 1000 if d_unit == "mm" else (d_val * 0.0254 if d_unit == "inch" else d_val)
    L = l_val * 1000 if l_unit == "km" else l_val
    area = math.pi * (D**2) / 4

    # [2. 유속 변환]
    v_calc_process = ""
    if v_unit == "m/s":
        v = v_val
    elif v_unit == "m³/s":
        v = v_val / area
    elif v_unit == "L/min":
        v = (v_val / 60000) / area
    elif v_unit == "BTU/hr (Water)":
        q_m3s = (v_val / 500) * 0.00006309
        v = q_m3s / area
        v_calc_process = "(BTU/hr 기반 환산)"

    # [3. Re 및 흐름 판정]
    re = (rho * v * D) / mu
    
    if re <= 2300:
        flow_type, f, color = "층류 (Laminar)", 64/re, "blue"
    elif 2300 < re < 4000:
        flow_type, f, color = "⚠️ 천이 영역 (Transition)", 0.316*(re**-0.25), "orange"
        st.warning(f"경고: Reynolds Number({re:.1f})가 천이 영역에 있습니다. 유동이 불안정합니다.")
    else:
        flow_type, f, color = "난류 (Turbulent)", 0.316*(re**-0.25), "red"

    # [4. 압력 손실 및 BTU 환산]
    sum_k = (n_elbow * 0.9) + (n_valve * 0.2)
    dp_major = f * (L / D) * (rho * v**2 / 2)
    dp_minor = sum_k * (rho * v**2 / 2)
    dp_total = dp_major + dp_minor
    
    energy_loss_btu_lb = (dp_total / rho) * 0.00042906

    # --- 결과 출력 ---
    st.divider()
    res_col1, res_col2 = st.columns(2)
    
    with res_col1:
        st.metric("Reynolds Number (Re)", f"{re:.1f}")
        st.metric("흐름 상태", flow_type)
        
    with res_col2:
        st.metric("총 압력 손실 (ΔP)", f"{dp_total/1000:.4f} kPa")
        st.metric("에너지 손실", f"{energy_loss_btu_lb:.6f} BTU/lb")

    with st.expander("📝 상세 계산 리포트 보기"):
        st.code(f"""
[상세 설계 리포트]
- 유체: {fluid_name} ({target_temp} °C)
- 관경(D): {D:.4f} m / 유속(v): {v:.4f} m/s
- 마찰계수(f): {f:.4f} ({flow_type})
- 주손실: {dp_major:.2f} Pa / 부손실: {dp_minor:.2f} Pa
- 최종 ΔP: {dp_total:.2f} Pa
        """)
