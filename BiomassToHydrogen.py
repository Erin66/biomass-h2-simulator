import pandas as pd
import streamlit as st

# ==============================================================================
# 核心计算逻辑 (与上一版相同，封装在类中)
# ==============================================================================
class BiomassToHydrogenCalculator:
    """
    一个根据输入原料和核心工艺参数，预测生物质制氢全流程产出的计算模型。
    """
    def __init__(self, rdf_feed_rate_daily, elemental_analysis, proximate_analysis, params):
        self.rdf_feed_rate = rdf_feed_rate_daily
        self.elemental = elemental_analysis
        self.proximate = proximate_analysis
        self.params = params
        
        self.atomic_weights = {'H': 1.008, 'C': 12.011, 'O': 15.999, 'N': 14.007, 'S': 32.06}
        self.molecular_weights = {
            'H2O': 18.015, 'CO': 28.01, 'H2': 2.016, 'CO2': 44.01,
            'CH4': 16.04, 'N2': 28.014, 'H2S': 34.08
        }
        self.results = {}

    def run_simulation(self):
        self._calculate_gasification()
        self._calculate_wgs()
        self._calculate_decarbonization()
        self._calculate_psa()
        self._calculate_steam_balance()
        return self.results

    def _calculate_gasification(self):
        oxygen_demand = self.rdf_feed_rate * self.params['oxygen_to_rdf_ratio']
        steam_demand = self.rdf_feed_rate * self.params['steam_to_rdf_ratio']
        
        C_in = self.rdf_feed_rate * (self.elemental['C'] / 100)
        H_in = self.rdf_feed_rate * (self.elemental['H'] / 100) + steam_demand * (2/18.015)
        O_in = self.rdf_feed_rate * (self.elemental['O'] / 100) + steam_demand * (16/18.015) + oxygen_demand
        N_in = self.rdf_feed_rate * (self.elemental['N'] / 100)
        S_in = self.rdf_feed_rate * (self.elemental['S'] / 100)
        Ash_in = self.rdf_feed_rate * (self.proximate['Ash'] / 100)

        unconverted_C = C_in * (1 - self.params['carbon_conversion_efficiency'])
        slag_output = Ash_in + unconverted_C
        
        gas_C = C_in - unconverted_C
        gas_H = H_in
        
        moles_C_total = gas_C / self.atomic_weights['C'] * 1000
        
        moles_CO = moles_C_total * self.params['syngas_mole_ratios']['CO']
        moles_CO2 = moles_C_total * self.params['syngas_mole_ratios']['CO2']
        moles_CH4 = moles_C_total * self.params['syngas_mole_ratios']['CH4']
        
        H_in_CH4 = moles_CH4 * 4 * self.atomic_weights['H'] / 1000
        moles_H2 = (gas_H - H_in_CH4) / self.molecular_weights['H2'] * 1000
        
        moles_N2 = N_in / self.molecular_weights['N2'] * 1000
        moles_H2S = S_in / self.molecular_weights['H2S'] * 1000
        
        syngas_composition_mass = {
            'H2': moles_H2 * self.molecular_weights['H2'] / 1000,
            'CO': moles_CO * self.molecular_weights['CO'] / 1000,
            'CO2': moles_CO2 * self.molecular_weights['CO2'] / 1000,
            'CH4': moles_CH4 * self.molecular_weights['CH4'] / 1000,
            'N2': moles_N2 * self.molecular_weights['N2'] / 1000,
            'H2S': moles_H2S * self.molecular_weights['H2S'] / 1000,
        }
        
        raw_syngas_mass = sum(syngas_composition_mass.values())
        
        self.results['gasification'] = {
            'oxygen_demand': oxygen_demand, 'steam_demand': steam_demand,
            'raw_syngas_mass': raw_syngas_mass, 'slag_output': slag_output,
            'syngas_composition_moles': {'H2': moles_H2, 'CO': moles_CO, 'CO2': moles_CO2, 
                                         'CH4': moles_CH4, 'N2': moles_N2, 'H2S': moles_H2S}
        }

    def _calculate_wgs(self):
        prev_moles = self.results['gasification']['syngas_composition_moles']
        co_to_convert = prev_moles['CO'] * self.params['wgs_co_conversion_rate']
        wgs_steam_demand_moles = co_to_convert * self.params['wgs_steam_to_co_ratio']
        wgs_steam_demand_mass = wgs_steam_demand_moles * self.molecular_weights['H2O'] / 1000
        
        shifted_moles = prev_moles.copy()
        shifted_moles['CO'] -= co_to_convert
        shifted_moles['H2'] += co_to_convert
        shifted_moles['CO2'] += co_to_convert
        
        shifted_syngas_mass = sum(v * self.molecular_weights.get(k, 0) for k, v in shifted_moles.items()) / 1000
        self.results['wgs'] = {'wgs_steam_demand': wgs_steam_demand_mass, 'shifted_syngas_mass': shifted_syngas_mass, 'shifted_syngas_moles': shifted_moles}

    def _calculate_decarbonization(self):
        prev_moles = self.results['wgs']['shifted_syngas_moles']
        co2_removed_moles = prev_moles['CO2'] * self.params['co2_removal_efficiency']
        co2_removed_mass = co2_removed_moles * self.molecular_weights['CO2'] / 1000
        
        decarbonized_moles = prev_moles.copy()
        decarbonized_moles['CO2'] -= co2_removed_moles
        decarbonized_syngas_mass = sum(v * self.molecular_weights.get(k, 0) for k, v in decarbonized_moles.items()) / 1000
        self.results['decarbonization'] = {'co2_product_mass': co2_removed_mass, 'decarbonized_syngas_mass': decarbonized_syngas_mass, 'decarbonized_syngas_moles': decarbonized_moles}

    def _calculate_psa(self):
        prev_moles = self.results['decarbonization']['decarbonized_syngas_moles']
        h2_in = prev_moles['H2']
        hydrogen_product_moles = h2_in * self.params['hydrogen_recovery_rate_psa']
        hydrogen_product_mass = hydrogen_product_moles * self.molecular_weights['H2'] / 1000
        psa_offgas_mass = self.results['decarbonization']['decarbonized_syngas_mass'] - hydrogen_product_mass
        self.results['psa'] = {'hydrogen_product_mass': hydrogen_product_mass, 'psa_offgas_mass': psa_offgas_mass}

    def _calculate_steam_balance(self):
        steam_gen_from_heat_recovery = self.rdf_feed_rate * self.params['steam_gen_per_ton_rdf']
        psa_offgas_lhv = 15
        heat_from_psa_offgas = self.results['psa']['psa_offgas_mass'] * 1000 * psa_offgas_lhv
        steam_gen_from_boiler = heat_from_psa_offgas / (2.5 * 1000)
        total_steam_generation = (steam_gen_from_heat_recovery + steam_gen_from_boiler) / 24
        
        steam_consum_mdea = self.results['decarbonization']['co2_product_mass'] * self.params['steam_per_ton_co2']
        steam_consum_wgs = self.results['wgs']['wgs_steam_demand']
        total_steam_consumption = (steam_consum_mdea + steam_consum_wgs) / 24
        self.results['steam_balance'] = {'total_steam_generation_t_h': total_steam_generation, 'total_steam_consumption_t_h': total_steam_consumption, 'net_steam_balance_t_h': total_steam_generation - total_steam_consumption}

# ==============================================================================
# Streamlit UI (用户交互界面)
# ==============================================================================

# 设置页面标题和图标
st.set_page_config(page_title="生物质制氢模拟器", page_icon="♻️")

# 标题
st.title("♻️ 生物质气化制氢全流程模拟器")
st.markdown("通过调整左侧的输入参数，实时预测制氢工厂的物料平衡和产出。")

# --- 侧边栏：所有用户输入控件 ---
st.sidebar.header("⚙️ 1. 输入原料参数")
rdf_feed = st.sidebar.number_input("每日RDF处理量 (t/d)", min_value=1.0, max_value=1000.0, value=60.0, step=1.0)

st.sidebar.subheader("元素分析 (%)")
c_perc = st.sidebar.number_input("碳 (C)", value=35.03)
h_perc = st.sidebar.number_input("氢 (H)", value=3.06)
o_perc = st.sidebar.number_input("氧 (O)", value=26.63)
n_perc = st.sidebar.number_input("氮 (N)", value=0.87)
s_perc = st.sidebar.number_input("硫 (S)", value=0.53)

st.sidebar.subheader("工业分析 (%)")
ash_perc = st.sidebar.number_input("灰分 (Ash)", value=33.57)


st.sidebar.header("🛠️ 2. 调整核心工艺参数")
# 使用滑块可以让用户更直观地感受参数变化带来的影响
st.sidebar.subheader("气化单元")
o2_rdf_ratio = st.sidebar.slider("氧气/RDF比 (t-O₂/t-RDF)", 0.1, 0.5, 0.23, 0.01)
steam_rdf_ratio = st.sidebar.slider("蒸汽/RDF比 (t-H₂O/t-RDF)", 0.1, 0.5, 0.228, 0.01)
c_conv = st.sidebar.slider("碳转化率 (%)", 90.0, 100.0, 99.9, 0.1) / 100.0

st.sidebar.subheader("变换与分离单元")
wgs_conv = st.sidebar.slider("CO变换转化率 (%)", 70.0, 99.0, 90.0, 1.0) / 100.0
h2_recovery = st.sidebar.slider("PSA氢气回收率 (%)", 70.0, 95.0, 85.0, 1.0) / 100.0

# --- 主页面：展示计算结果 ---
st.header("📊 计算结果总览")

# 创建一个按钮来触发计算
if st.button("运行模拟计算"):
    # 收集UI界面的输入数据
    user_elemental = {'C': c_perc, 'H': h_perc, 'O': o_perc, 'N': n_perc, 'S': s_perc}
    user_proximate = {'Ash': ash_perc}
    user_params = {
        'oxygen_to_rdf_ratio': o2_rdf_ratio,
        'steam_to_rdf_ratio': steam_rdf_ratio,
        'carbon_conversion_efficiency': c_conv,
        'syngas_mole_ratios': {'CO': 0.25, 'CO2': 0.35, 'CH4': 0.05}, # 简化模型，保持固定
        'wgs_co_conversion_rate': wgs_conv,
        'wgs_steam_to_co_ratio': 1.5,
        'co2_removal_efficiency': 0.98,
        'hydrogen_recovery_rate_psa': h2_recovery,
        'steam_gen_per_ton_rdf': 1.27,
        'steam_per_ton_co2': 0.5
    }

    # 运行计算
    with st.spinner('正在进行全流程模拟计算...'):
        calculator = BiomassToHydrogenCalculator(
            rdf_feed_rate_daily=rdf_feed,
            elemental_analysis=user_elemental,
            proximate_analysis=user_proximate,
            params=user_params
        )
        results = calculator.run_simulation()
    
    st.success('计算完成！')

    # 使用列布局来美化结果展示
    col1, col2 = st.columns(2)

    # **核心产品**
    with col1:
        st.subheader("⭐ 核心产品产量")
        st.metric(label="最终氢气产量 (t/d)", value=f"{results['psa']['hydrogen_product_mass']:.3f}")
    
    # **主要输入**
    with col2:
        st.subheader("💧 主要输入需求")
        st.metric(label="总蒸汽需求 (t/d)", value=f"{results['gasification']['steam_demand'] + results['wgs']['wgs_steam_demand']:.2f}")
        st.metric(label="纯氧需求 (t/d)", value=f"{results['gasification']['oxygen_demand']:.2f}")

    st.markdown("---")

    # **详细物料平衡**
    st.subheader("📑 详细物料平衡 (t/d)")
    
    # 创建一个数据框来展示
    balance_data = {
        "物料流": ["粗合成气", "变换后合成气", "富氢合成气 (入PSA)", "PSA尾气 (燃料)", "气化灰渣", "CO₂产品"],
        "产量 (t/d)": [
            results['gasification']['raw_syngas_mass'],
            results['wgs']['shifted_syngas_mass'],
            results['decarbonization']['decarbonized_syngas_mass'],
            results['psa']['psa_offgas_mass'],
            results['gasification']['slag_output'],
            results['decarbonization']['co2_product_mass']
        ]
    }
    df_balance = pd.DataFrame(balance_data).round(2)
    st.dataframe(df_balance)

    st.markdown("---")

    # **蒸汽能量平衡**
    st.subheader("♨️ 全厂蒸汽能量平衡")
    steam_gen = results['steam_balance']['total_steam_generation_t_h']
    steam_con = results['steam_balance']['total_steam_consumption_t_h']
    steam_net = results['steam_balance']['net_steam_balance_t_h']
    
    col_s1, col_s2, col_s3 = st.columns(3)
    col_s1.metric("蒸汽总产量 (t/h)", f"{steam_gen:.2f}")
    col_s2.metric("蒸汽总消耗 (t/h)", f"{steam_con:.2f}")
    col_s3.metric("净蒸汽平衡 (t/h)", f"{steam_net:.2f}", delta=f"{steam_net:.2f}")

    # 可视化图表
    steam_df = pd.DataFrame({'类型': ['产量', '消耗'], '蒸汽量 (t/h)': [steam_gen, steam_con]})
    st.bar_chart(steam_df.set_index('类型'))
    
    if steam_net >= 0:
        st.info("💡 结论: 系统能量可自给自足，并有富余蒸汽。")
    else:
        st.warning("⚠️ 结论: 系统能量存在缺口，需要外部补充蒸汽。")

else:
    st.info("⬅️ 请在左侧配置参数，然后点击“运行模拟计算”按钮。")