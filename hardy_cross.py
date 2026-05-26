import streamlit as st
import math
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd

# --- 1. 파이프 객체 정의 ---
class Pipe:
    def __init__(self, id, node_start, node_end, length, diameter, roughness, initial_q=0.0):
        self.id = id
        self.start = node_start
        self.end = node_end
        self.L = length
        self.D = diameter
        self.epsilon = roughness
        self.Q = initial_q  
        self.rho = 1000
        self.mu = 0.89e-3
        
    def get_c_constant(self):
        return (8 * self.rho * self.L) / (math.pi**2 * (self.D**5))

    def get_friction_factor(self):
        re = (4 * self.rho * abs(self.Q)) / (math.pi * self.D * self.mu)
        if re < 2300: 
            return 64 / re if re > 0 else 0.02
        term1 = self.epsilon / (3.7065 * self.D)
        term2 = 5.0452 / re
        term3 = (1 / 2.8257) * (self.epsilon / self.D)**1.1098 + 5.8506 / (re**0.8981)
        return (-2.0 * math.log10(term1 - term2 * math.log10(term3)))**-2

    def get_delta_p(self):
        return self.get_c_constant() * self.get_friction_factor() * self.Q * abs(self.Q)

    def get_dp_over_q(self):
        return abs(self.get_delta_p() / self.Q) if self.Q != 0 else 0

# --- 2. 동적 하디 크로스 수치 해석 엔진 ---
def run_dynamic_hardy_cross(pipes, loops_nodes, max_iter=30, tolerance=1e-5):
    iteration_history = []
    pipe_dict = {(p.start, p.end): p for p in pipes}
    pipe_dict.update({(p.end, p.start): p for p in pipes})
    
    for idx in range(max_iter):
        dq_list = []
        
        for loop in loops_nodes:
            sum_dp = 0.0
            sum_dp_q = 0.0
            
            for i in range(len(loop)):
                n1 = loop[i]
                n2 = loop[(i + 1) % len(loop)]
                
                if (n1, n2) in pipe_dict:
                    p = pipe_dict[(n1, n2)]
                    sign = 1.0 if p.start == n1 else -1.0
                    sum_dp += sign * p.get_delta_p()
                    sum_dp_q += p.get_dp_over_q()
            
            dq = -sum_dp / (2 * sum_dp_q) if sum_dp_q != 0 else 0
            dq_list.append(dq)
            
        q_adjustments = {p.id: 0.0 for p in pipes}
        for l_idx, loop in enumerate(loops_nodes):
            dq = dq_list[l_idx]
            for i in range(len(loop)):
                n1 = loop[i]
                n2 = loop[(i + 1) % len(loop)]
                p = pipe_dict[(n1, n2)]
                sign = 1.0 if p.start == n1 else -1.0
                q_adjustments[p.id] += sign * dq
                
        for p in pipes:
            p.Q += q_adjustments[p.id]
            
        iteration_history.append(dq_list)
        if all(abs(dq) < tolerance for dq in dq_list):
            break
            
    return iteration_history

# --- 3. Streamlit 실행 메인 함수 ---
def run_hardy_cross():
    if 'pipe_data' not in st.session_state:
        st.session_state.pipe_data = [
            {"id": 1, "start": "A", "end": "B", "L": 300.0, "D": 0.250, "init_q": 0.060, "sx": 2.0, "sy": 0.0, "ex": 1.0, "ey": 0.0},
            {"id": 2, "start": "B", "end": "E", "L": 250.0, "D": 0.200, "init_q": 0.020, "sx": 1.0, "sy": 0.0, "ex": 1.0, "ey": 1.0},
            {"id": 3, "start": "E", "end": "D", "L": 300.0, "D": 0.200, "init_q": -0.040, "sx": 1.0, "sy": 1.0, "ex": 2.0, "ey": 1.0},
            {"id": 4, "start": "D", "end": "A", "L": 250.0, "D": 0.250, "init_q": -0.065, "sx": 2.0, "sy": 1.0, "ex": 2.0, "ey": 0.0},
            {"id": 5, "start": "B", "end": "C", "L": 300.0, "D": 0.200, "init_q": 0.040, "sx": 1.0, "sy": 0.0, "ex": 0.0, "ey": 0.0},
            {"id": 6, "start": "C", "end": "F", "L": 250.0, "D": 0.200, "init_q": 0.028, "sx": 0.0, "sy": 0.0, "ex": 0.0, "ey": 1.0},
            {"id": 7, "start": "F", "end": "E", "L": 300.0, "D": 0.150, "init_q": -0.035, "sx": 0.0, "sy": 1.0, "ex": 1.0, "ey": 1.0}
        ]

    st.sidebar.header("[1] 시스템 가동 조건 설정")
    control_node = st.sidebar.selectbox("펌프/밸브 제어 노드 선택", ["A", "B", "C", "D", "E", "F"])
    total_inflow = st.sidebar.slider("시스템 총 유입 유량 (m³/s)", 0.05, 0.40, 0.125, step=0.005)
    roughness_val = st.sidebar.number_input("관 절대 조도 (m)", value=0.00025, format="%.5f")
    pump_eff = st.sidebar.slider("펌프 효율 (η)", 0.5, 0.9, 0.75)

    st.subheader("➕ 설계 배관망 동적 추가(절대 좌표)")
    
    with st.expander("📐 새로운 파이프라인 및 XY 좌표 지정하여 추가하기", expanded=False):
        c1, c2, c3, c4, c5 = st.columns(5)
        new_id = len(st.session_state.pipe_data) + 1
        p_start = c1.text_input("시작 노드 (예: G)", value="G").strip().upper()
        p_end = c2.text_input("끝 노드 (예: H)", value="H").strip().upper()
        p_L = c3.number_input("길이 L (m)", min_value=1.0, value=200.0)
        p_D = c4.number_input("직경 D (m)", min_value=0.01, value=0.200, format="%.3f")
        p_q = c5.number_input("가정 유량", value=0.010, format="%.3f")
        
        cx1, cy1, cx2, cy2 = st.columns(4)
        sx = cx1.number_input("시작 노드 X 위치", value=3.0, step=1.0)
        sy = cy1.number_input("시작 노드 Y 위치", value=0.0, step=1.0)
        ex = cx2.number_input("끝 노드 X 위치", value=3.0, step=1.0)
        ey = cy2.number_input("끝 노드 Y 위치", value=1.0, step=1.0)
        
        if st.button("➕ 네트워크에 배관 및 좌표 등록", use_container_width=True):
            if p_start == p_end:
                st.error("시작 노드와 끝 노드는 같을 수 없습니다.")
            else:
                st.session_state.pipe_data.append({
                    "id": new_id, "start": p_start, "end": p_end, "L": p_L, "D": p_D, "init_q": p_q,
                    "sx": sx, "sy": sy, "ex": ex, "ey": ey
                })
                st.success(f"Pipe {new_id} ({p_start}➔{p_end})가 완벽히 배치되었습니다!")
                st.rerun()

    df_pipes = pd.DataFrame(st.session_state.pipe_data)
    st.dataframe(df_pipes, use_container_width=True)
    

    reset_col1, reset_col2 = st.columns(2)
    
    if reset_col1.button("🔄 교재 예제 데이터로 초기화", use_container_width=True):
        # 세션 데이터를 완전히 지운 후 새로고침하여 초기 예제 배열(그림 5.2) 복원
        if 'pipe_data' in st.session_state:
            del st.session_state.pipe_data
        st.success("교재 표준 그리드 배관망 데이터가 복원되었습니다.")
        st.rerun()
        
    if reset_col2.button("🗑️ 전체 노드 삭제 (Blank Reset)", use_container_width=True):
        # 모든 배관 제거
        st.session_state.pipe_data = []
        st.warning("네트워크의 모든 노드와 배관이 삭제되었습니다. 새로운 설계를 시작하세요.")
        st.rerun()

    G_setup = nx.Graph()
    pipes = []
    pos = {}  
    
    scale = total_inflow / 0.125
    for p in st.session_state.pipe_data:
        pipes.append(Pipe(p['id'], p['start'], p['end'], p['L'], p['D'], roughness_val, p['init_q'] * scale))
        G_setup.add_edge(p['start'], p['end'])
        
        pos[p['start']] = (p.get('sx', 0.0), p.get('sy', 0.0))
        pos[p['end']] = (p.get('ex', 0.0), p.get('ey', 0.0))

    auto_loops = nx.cycle_basis(G_setup)

    st.info(f"🔍 시스템 위상 분석 완료: 감지된 독립 루프 개수 = **{len(auto_loops)}개**")

    if len(auto_loops) == 0:
        st.warning("⚠️ 현재 폐회로(Loop)가 없는 개방형 계통입니다. 폐회로가 구성되도록 노드를 이어주세요.")
        return

    history = run_dynamic_hardy_cross(pipes, auto_loops)

    col1, col2 = st.columns([3, 2])
    
    with col1:
        st.subheader("🖼️ 절대 좌표 정렬 기반 2D 배관 플랜 도면")
        
        G_draw = nx.DiGraph()
        for p in pipes:
            if p.Q >= 0: G_draw.add_edge(p.start, p.end, weight=abs(p.Q), id=p.id)
            else: G_draw.add_edge(p.end, p.start, weight=abs(p.Q), id=p.id)
            
        fig, ax = plt.subplots(figsize=(7, 4.5))
        
        nx.draw_networkx_nodes(G_draw, pos, node_size=600, node_color='#D6EAF8', ax=ax)
        nx.draw_networkx_labels(G_draw, pos, font_size=11, font_weight='bold', ax=ax)
        
        if control_node in pos:
            nx.draw_networkx_nodes(G_draw, pos, nodelist=[control_node], node_size=700, node_color='#FF5733', ax=ax)
            
        edges = G_draw.edges(data=True)
        weights = [max(e[2]['weight'] * 150, 1.5) for e in edges]
        nx.draw_networkx_edges(G_draw, pos, width=weights, edge_color='#2C3E50', arrowsize=18, ax=ax)
        
        edge_labels = {(p.start, p.end) if p.Q >= 0 else (p.end, p.start): f"#{p.id}\n{abs(p.Q):.3f}m³/s" for p in pipes}
        nx.draw_networkx_edge_labels(G_draw, pos, edge_labels=edge_labels, font_size=8, font_color='red')
        
        plt.tight_layout()
        st.pyplot(fig)

    with col2:
        st.subheader("📊 자동 다중루프 보정 수렴 리포트")
        
        report_data = []
        for idx, dqs in enumerate(history[:5]):
            row = {"반복 횟수": f"{idx+1}차"}
            for l_idx, dq_val in enumerate(dqs):
                row[f"Loop {l_idx+1} ΔQ"] = f"{dq_val:.6f}"
            report_data.append(row)
        st.table(report_data)
        
        total_dp_loss = sum(abs(p.get_delta_p()) for p in pipes)
        required_power_w = total_dp_loss * total_inflow
        required_power_kw = required_power_w / (pump_eff * 1000)
        
        st.metric("시스템 총 마찰 압력 강하", f"{total_dp_loss:,.1f} N/m²")
        st.metric("💡 권장 최소 펌프 정격 동력", f"{required_power_kw:.2f} kW")

    st.divider()
    
    st.subheader("📋 파이프 라인별 해석 결과 상세 내역")
    result_table = []
    for p in pipes:
        result_table.append({
            "배관 번호": f"Pipe {p.id}",
            "연결 노드": f"{p.start} ➔ {p.end}",
            "최종 유량 (m³/s)": f"{p.Q:.4f}",
            "마찰 계수 (f)": f"{p.get_friction_factor():.4f}",
            "압력 강하 (N/m²)": f"{p.get_delta_p():,.1f}"
        })
    st.dataframe(result_table, use_container_width=True)
