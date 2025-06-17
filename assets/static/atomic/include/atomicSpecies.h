#ifndef SPECIES_MASS_DATA_H
#define SPECIES_MASS_DATA_H
#include <unordered_map>
#include <string_view>
#include <string>

namespace chemSpecies {
    struct Species {
        const std::string_view m_name; //< Name of the species
        const std::string_view m_el; //< Element symbol
        const int m_nz; //< NZ
        const int m_n; //< N
        const int m_z; //< Z
        const int m_a; //< A
        const double m_bindingEnergy; //< Binding energy
        const std::string_view m_betaCode; //< Beta decay code
        const double m_betaDecayEnergy; //< Beta decay energy
        const double m_atomicMass; //< Atomic mass
        const double m_atomicMassUnc; //< Atomic mass uncertainty

        Species(const std::string_view name, const std::string_view el, const int nz, const int n, const int z, const int a, const double bindingEnergy, const std::string_view betaCode, const double betaDecayEnergy, const double atomicMass, const double atomicMassUnc)
            : m_name(name), m_el(el), m_nz(nz), m_n(n), m_z(z), m_a(a), m_bindingEnergy(bindingEnergy), m_betaCode(betaCode), m_betaDecayEnergy(betaDecayEnergy), m_atomicMass(atomicMass), m_atomicMassUnc(atomicMassUnc) {};

        Species(const Species& species) {
            m_name = species.m_name;
            m_el = species.m_el;
            m_nz = species.m_nz;
            m_n = species.m_n;
            m_z = species.m_z;
            m_a = species.m_a;
            m_bindingEnergy = species.m_bindingEnergy;
            m_betaCode = species.m_betaCode;
            m_betaDecayEnergy = species.m_betaDecayEnergy;
            m_atomicMass = species.m_atomicMass;
            m_atomicMassUnc = species.m_atomicMassUnc;
        }

        double mass() const {
            return m_atomicMass;
        }

        double massUnc() const {
            return m_atomicMassUnc;
        }

        double bindingEnergy() const {
            return m_bindingEnergy;
        }

        double betaDecayEnergy() const {
            return m_betaDecayEnergy;
        }

        std::string_view betaCode() const {
            return m_betaCode;
        }

        std::string_view name() const {
            return m_name;
        }

        std::string_view el() const {
            return m_el;
        }

        int nz() const {
            return m_nz;
        }

        int n() const {    
            return m_n;
        }

        int z() const {
            return m_z;
        }

        int a() const {
            return m_a;
        }

        friend std::ostream& operator<<(std::ostream& os, const Species& species) {
            os << static_cast<std::string>(species.m_name) << " (" << species.m_atomicMass << " u)";
            return os;
        }
    };
    static const Species n_1("n-1", "n", 1, 1, 0, 1, 0.0, "B-", 782.347, 1.0086649159, 0.00047);
    static const Species H_1("H-1", "H", -1, 0, 1, 1, 0.0, "B-", std::numeric_limits<double>::quiet_NaN(), 1.007825031898, 1.4e-05);
    static const Species H_2("H-2", "H", 0, 1, 1, 2, 1112.2831, "B-", std::numeric_limits<double>::quiet_NaN(), 2.014101777844, 1.5e-05);
    static const Species H_3("H-3", "H", 1, 2, 1, 3, 2827.2654, "B-", 18.59202, 3.01604928132, 8e-05);
    static const Species He_3("He-3", "He", -1, 1, 2, 3, 2572.68044, "B-", -13736.0, 3.01602932197, 6e-05);
    static const Species Li_3("Li-3", "Li", -3, 0, 3, 3, -2267.0, "B-", std::numeric_limits<double>::quiet_NaN(), 3.030775, 2147.0);
    static const Species H_4("H-4", "H", 2, 3, 1, 4, 1720.4491, "B-", 22196.2131, 4.026431867, 107.354);
    static const Species He_4("He-4", "He", 0, 2, 2, 4, 7073.9156, "B-", -22898.274, 4.00260325413, 0.00016);
    static const Species Li_4("Li-4", "Li", -2, 1, 3, 4, 1153.7603, "B-", std::numeric_limits<double>::quiet_NaN(), 4.027185561, 227.733);
    static const Species H_5("H-5", "H", 3, 4, 1, 5, 1336.3592, "B-", 21661.2131, 5.035311492, 96.02);
    static const Species He_5("He-5", "He", 1, 3, 2, 5, 5512.1325, "B-", -447.6529, 5.012057224, 21.47);
    static const Species Li_5("Li-5", "Li", -1, 2, 3, 5, 5266.1325, "B-", -25460.0, 5.0125378, 53.677);
    static const Species Be_5("Be-5", "Be", -3, 1, 4, 5, 18.0, "B-", std::numeric_limits<double>::quiet_NaN(), 5.03987, 2150.0);
    static const Species H_6("H-6", "H", 4, 5, 1, 6, 961.6395, "B-", 24283.6294, 6.044955437, 272.816);
    static const Species He_6("He-6", "He", 2, 4, 2, 6, 4878.5199, "B-", 3505.2147, 6.018885889, 0.057);
    static const Species Li_6("Li-6", "Li", 0, 3, 3, 6, 5332.3312, "B-", -4288.1534, 6.01512288742, 0.00155);
    static const Species Be_6("Be-6", "Be", -2, 2, 4, 6, 4487.2478, "B-", -28945.0, 6.0197264090000004, 5.848);
    static const Species B_6("B-6", "B", -4, 1, 5, 6, -467.0, "B-", std::numeric_limits<double>::quiet_NaN(), 6.0508, 2150.0);
    static const Species H_7("H-7", "H", 5, 6, 1, 7, 940.0, "B-", 23062.0, 7.052749, 1078.0);
    static const Species He_7("He-7", "He", 3, 5, 2, 7, 4123.0578, "B-", 11166.0229, 7.027990652, 8.115);
    static const Species Li_7("Li-7", "Li", 1, 4, 3, 7, 5606.4401, "B-", -861.893, 7.01600343426, 0.0045);
    static const Species Be_7("Be-7", "Be", -1, 3, 4, 7, 5371.5487, "B-", -11907.5551, 7.016928714, 0.076);
    static const Species B_7("B-7", "B", -3, 2, 5, 7, 3558.7055, "B-", std::numeric_limits<double>::quiet_NaN(), 7.029712, 27.0);
    static const Species He_8("He-8", "He", 4, 6, 2, 8, 3924.521, "B-", 10663.8784, 8.033934388, 0.095);
    static const Species Li_8("Li-8", "Li", 2, 5, 3, 8, 5159.7124, "B-", 16004.1329, 8.022486244, 0.05);
    static const Species Be_8("Be-8", "Be", 0, 4, 4, 8, 7062.4356, "B-", -17979.8973, 8.005305102, 0.037);
    static const Species B_8("B-8", "B", -2, 3, 5, 8, 4717.1551, "B-", -12142.7002, 8.024607315, 1.073);
    static const Species C_8("C-8", "C", -4, 2, 6, 8, 3101.5242, "B-", std::numeric_limits<double>::quiet_NaN(), 8.037643039, 19.584);
    static const Species He_9("He-9", "He", 5, 7, 2, 9, 3349.038, "B-", 15980.9213, 9.043946414, 50.259);
    static const Species Li_9("Li-9", "Li", 3, 6, 3, 9, 5037.7685, "B-", 13606.4541, 9.026790191, 0.2);
    static const Species Be_9("Be-9", "Be", 1, 5, 4, 9, 6462.6693, "B-", -1068.0349, 9.012183062, 0.082);
    static const Species B_9("B-9", "B", -1, 4, 5, 9, 6257.0713, "B-", -16494.4854, 9.013329645, 0.969);
    static const Species C_9("C-9", "C", -3, 3, 6, 9, 4337.4233, "B-", std::numeric_limits<double>::quiet_NaN(), 9.031037202, 2.293);
    static const Species He_10("He-10", "He", 6, 8, 2, 10, 2995.134, "B-", 16144.5191, 10.052815306, 99.676);
    static const Species Li_10("Li-10", "Li", 4, 7, 3, 10, 4531.3512, "B-", 20445.1411, 10.035483453, 13.656);
    static const Species Be_10("Be-10", "Be", 2, 6, 4, 10, 6497.6306, "B-", 556.8759, 10.013534692, 0.086);
    static const Species B_10("B-10", "B", 0, 5, 5, 10, 6475.0835, "B-", -3648.0623, 10.012936862, 0.016);
    static const Species C_10("C-10", "C", -2, 4, 6, 10, 6032.0426, "B-", -23101.3545, 10.016853217, 0.075);
    static const Species N_10("N-10", "N", -4, 3, 7, 10, 3643.6724, "B-", std::numeric_limits<double>::quiet_NaN(), 10.04165354, 429.417);
    static const Species Li_11("Li-11", "Li", 5, 8, 3, 11, 4155.3817, "B-", 20551.0898, 11.043723581, 0.66);
    static const Species Be_11("Be-11", "Be", 3, 7, 4, 11, 5952.5402, "B-", 11509.4607, 11.02166108, 0.255);
    static const Species B_11("B-11", "B", 1, 6, 5, 11, 6927.7323, "B-", -1981.6889, 11.009305166, 0.013);
    static const Species C_11("C-11", "C", -1, 5, 6, 11, 6676.4563, "B-", -13716.2469, 11.011432597, 0.064);
    static const Species N_11("N-11", "N", -3, 4, 7, 11, 5358.4023, "B-", -23373.2693, 11.026157593, 5.368);
    static const Species O_11("O-11", "O", -5, 3, 8, 11, 3162.4372, "B-", std::numeric_limits<double>::quiet_NaN(), 11.051249828, 64.453);
    static const Species Li_12("Li-12", "Li", 6, 9, 3, 12, 3791.5999, "B-", 23931.8152, 12.052613942, 32.213);
    static const Species Be_12("Be-12", "Be", 4, 8, 4, 12, 5720.7223, "B-", 11708.3636, 12.026922082, 2.048);
    static const Species B_12("B-12", "B", 2, 7, 5, 12, 6631.2237, "B-", 13369.3979, 12.014352638, 1.418);
    static const Species C_12("C-12", "C", 0, 6, 6, 12, 7680.1446, "B-", -17338.0681, 12.0, 0.0);
    static const Species N_12("N-12", "N", -2, 5, 7, 12, 6170.11, "B-", -14675.2668, 12.01861318, 1.073);
    static const Species O_12("O-12", "O", -4, 4, 8, 12, 4881.9755, "B-", std::numeric_limits<double>::quiet_NaN(), 12.034367726, 12.882);
    static const Species Li_13("Li-13", "Li", 7, 10, 3, 13, 3507.6307, "B-", 23321.8152, 13.061171503, 75.15);
    static const Species Be_13("Be-13", "Be", 5, 9, 4, 13, 5241.4359, "B-", 17097.1315, 13.036134506, 10.929);
    static const Species B_13("B-13", "B", 3, 8, 5, 13, 6496.4194, "B-", 13436.9387, 13.017779981, 1.073);
    static const Species C_13("C-13", "C", 1, 7, 6, 13, 7469.8495, "B-", -2220.4718, 13.00335483534, 0.00025);
    static const Species N_13("N-13", "N", -1, 6, 7, 13, 7238.8634, "B-", -17769.9506, 13.005738609, 0.289);
    static const Species O_13("O-13", "O", -3, 5, 8, 13, 5811.7636, "B-", -18915.0, 13.024815435, 10.226);
    static const Species F_13("F-13", "F", -5, 4, 9, 13, 4297.0, "B-", std::numeric_limits<double>::quiet_NaN(), 13.045121, 537.0);
    static const Species Be_14("Be-14", "Be", 6, 10, 4, 14, 4993.8973, "B-", 16290.8166, 14.04289292, 141.97);
    static const Species B_14("B-14", "B", 4, 9, 5, 14, 6101.6451, "B-", 20643.7926, 14.02540401, 22.773);
    static const Species C_14("C-14", "C", 2, 8, 6, 14, 7520.3198, "B-", 156.4765, 14.00324198862, 0.00403);
    static const Species N_14("N-14", "N", 0, 7, 7, 14, 7475.6148, "B-", -5144.3643, 14.00307400425, 0.00024);
    static const Species O_14("O-14", "O", -2, 6, 8, 14, 7052.2783, "B-", -23956.6215, 14.008596706, 0.027);
    static const Species F_14("F-14", "F", -4, 5, 9, 14, 5285.2091, "B-", std::numeric_limits<double>::quiet_NaN(), 14.034315196, 44.142);
    static const Species Be_15("Be-15", "Be", 7, 11, 4, 15, 4540.9708, "B-", 20868.4411, 15.053490215, 177.99);
    static const Species B_15("B-15", "B", 5, 10, 5, 15, 5880.0438, "B-", 19084.2343, 15.031087023, 22.575);
    static const Species C_15("C-15", "C", 3, 9, 6, 15, 7100.1696, "B-", 9771.7071, 15.010599256, 0.858);
    static const Species N_15("N-15", "N", 1, 8, 7, 15, 7699.4603, "B-", -2754.1841, 15.00010889827, 0.00062);
    static const Species O_15("O-15", "O", -1, 7, 8, 15, 7463.6915, "B-", -13711.13, 15.003065636, 0.526);
    static const Species F_15("F-15", "F", -3, 6, 9, 15, 6497.4597, "B-", -23648.6215, 15.017785139, 15.029);
    static const Species Ne_15("Ne-15", "Ne", -5, 5, 10, 15, 4868.7285, "B-", std::numeric_limits<double>::quiet_NaN(), 15.043172977, 71.588);
    static const Species Be_16("Be-16", "Be", 8, 12, 4, 16, 4285.2851, "B-", 20335.4399, 16.061672036, 177.99);
    static const Species B_16("B-16", "B", 6, 11, 5, 16, 5507.3535, "B-", 23417.5656, 16.039841045, 26.373);
    static const Species C_16("C-16", "C", 4, 10, 6, 16, 6922.0546, "B-", 8010.226, 16.014701255, 3.84);
    static const Species N_16("N-16", "N", 2, 9, 7, 16, 7373.7971, "B-", 10420.9094, 16.006101925, 2.47);
    static const Species O_16("O-16", "O", 0, 8, 8, 16, 7976.2072, "B-", -15412.184, 15.99491461926, 0.00032);
    static const Species F_16("F-16", "F", -2, 7, 9, 16, 6964.049, "B-", -13311.5932, 16.011460278, 5.758);
    static const Species Ne_16("Ne-16", "Ne", -4, 6, 10, 16, 6083.1777, "B-", std::numeric_limits<double>::quiet_NaN(), 16.02575086, 21.986);
    static const Species B_17("B-17", "B", 7, 12, 5, 17, 5269.6677, "B-", 22684.4422, 17.046931399, 219.114);
    static const Species C_17("C-17", "C", 5, 11, 6, 17, 6558.0262, "B-", 13161.8007, 17.02257865, 18.641);
    static const Species N_17("N-17", "N", 3, 10, 7, 17, 7286.2294, "B-", 8678.843, 17.008448876, 16.103);
    static const Species O_17("O-17", "O", 1, 9, 8, 17, 7750.7291, "B-", -2760.4655, 16.99913175595, 0.00069);
    static const Species F_17("F-17", "F", -1, 8, 9, 17, 7542.3284, "B-", -14548.7507, 17.002095237, 0.266);
    static const Species Ne_17("Ne-17", "Ne", -3, 7, 10, 17, 6640.4991, "B-", -18219.1277, 17.017713962, 0.38);
    static const Species Na_17("Na-17", "Na", -5, 6, 11, 17, 5522.7653, "B-", std::numeric_limits<double>::quiet_NaN(), 17.037273, 64.0);
    static const Species B_18("B-18", "B", 8, 13, 5, 18, 4976.6306, "B-", 26873.3742, 18.055601683, 219.18);
    static const Species C_18("C-18", "C", 6, 12, 6, 18, 6426.1321, "B-", 11806.0982, 18.02675193, 32.206);
    static const Species N_18("N-18", "N", 4, 11, 7, 18, 7038.5627, "B-", 13895.9838, 18.014077563, 19.935);
    static const Species O_18("O-18", "O", 2, 10, 8, 18, 7767.0981, "B-", -1655.9288, 17.99915961214, 0.00069);
    static const Species F_18("F-18", "F", 0, 9, 9, 18, 7631.6383, "B-", -4444.5049, 18.000937324, 0.497);
    static const Species Ne_18("Ne-18", "Ne", -2, 8, 10, 18, 7341.2577, "B-", -19720.3745, 18.005708696, 0.39);
    static const Species Na_18("Na-18", "Na", -4, 7, 11, 18, 6202.2176, "B-", std::numeric_limits<double>::quiet_NaN(), 18.026879388, 100.785);
    static const Species B_19("B-19", "B", 9, 14, 5, 19, 4719.6346, "B-", 27356.4961, 19.064166, 564.0);
    static const Species C_19("C-19", "C", 7, 13, 6, 19, 6118.274, "B-", 16557.4995, 19.034797594, 105.625);
    static const Species N_19("N-19", "N", 5, 12, 7, 19, 6948.5452, "B-", 12523.3972, 19.017022389, 17.61);
    static const Species O_19("O-19", "O", 3, 11, 8, 19, 7566.4952, "B-", 4820.3029, 19.003577969, 2.83);
    static const Species F_19("F-19", "F", 1, 10, 9, 19, 7779.0192, "B-", -3239.4986, 18.99840316207, 0.00088);
    static const Species Ne_19("Ne-19", "Ne", -1, 9, 10, 19, 7567.3431, "B-", -11177.331, 19.001880906, 0.171);
    static const Species Na_19("Na-19", "Na", -3, 8, 11, 19, 6937.8864, "B-", -18909.0095, 19.013880264, 11.309);
    static const Species Mg_19("Mg-19", "Mg", -5, 7, 12, 19, 5901.4992, "B-", std::numeric_limits<double>::quiet_NaN(), 19.03417992, 64.413);
    static const Species B_20("B-20", "B", 10, 15, 5, 20, 4405.6529, "B-", 31898.0019, 20.074505644, 586.538);
    static const Species C_20("C-20", "C", 8, 14, 6, 20, 5961.4356, "B-", 15737.0689, 20.040261732, 247.585);
    static const Species N_20("N-20", "N", 6, 13, 7, 20, 6709.1717, "B-", 17970.3261, 20.023367295, 84.696);
    static const Species O_20("O-20", "O", 4, 12, 8, 20, 7568.5707, "B-", 3813.6349, 20.004075357, 0.95);
    static const Species F_20("F-20", "F", 2, 11, 9, 20, 7720.1351, "B-", 7024.4689, 19.999981252, 0.031);
    static const Species Ne_20("Ne-20", "Ne", 0, 10, 10, 20, 8032.2412, "B-", -13892.4207, 19.99244017525, 0.00165);
    static const Species Na_20("Na-20", "Na", -2, 9, 11, 20, 7298.5028, "B-", -10627.2054, 20.007354301, 1.19);
    static const Species Mg_20("Mg-20", "Mg", -4, 8, 12, 20, 6728.0252, "B-", std::numeric_limits<double>::quiet_NaN(), 20.018763075, 2.0);
    static const Species B_21("B-21", "B", 11, 16, 5, 21, 4152.5265, "B-", 32740.0, 21.084147485, 599.75);
    static const Species C_21("C-21", "C", 9, 15, 6, 21, 5674.0, "B-", 20411.0, 21.049, 640.0);
    static const Species N_21("N-21", "N", 7, 14, 7, 21, 6609.0159, "B-", 17169.8816, 21.027087573, 143.906);
    static const Species O_21("O-21", "O", 5, 13, 8, 21, 7389.3747, "B-", 8109.639, 21.008654948, 12.882);
    static const Species F_21("F-21", "F", 3, 12, 9, 21, 7738.2934, "B-", 5684.1712, 20.999948893, 1.932);
    static const Species Ne_21("Ne-21", "Ne", 1, 11, 10, 21, 7971.7136, "B-", -3546.919, 20.993846685, 0.041);
    static const Species Na_21("Na-21", "Na", -1, 10, 11, 21, 7765.5581, "B-", -13088.708, 20.997654459, 0.045);
    static const Species Mg_21("Mg-21", "Mg", -3, 9, 12, 21, 7105.0317, "B-", -16186.0, 21.011705764, 0.81);
    static const Species Al_21("Al-21", "Al", -5, 8, 13, 21, 6297.0, "B-", std::numeric_limits<double>::quiet_NaN(), 21.029082, 644.0);
    static const Species C_22("C-22", "C", 10, 16, 6, 22, 5421.0778, "B-", 21846.3983, 22.05755399, 248.515);
    static const Species N_22("N-22", "N", 8, 15, 7, 22, 6378.5347, "B-", 22481.7725, 22.034100918, 223.06);
    static const Species O_22("O-22", "O", 6, 14, 8, 22, 7364.8722, "B-", 6489.6562, 22.009965744, 61.107);
    static const Species F_22("F-22", "F", 4, 13, 9, 22, 7624.2954, "B-", 10818.0916, 22.002998812, 13.31);
    static const Species Ne_22("Ne-22", "Ne", 2, 12, 10, 22, 8080.4656, "B-", -2843.3243, 21.991385113, 0.018);
    static const Species Na_22("Na-22", "Na", 0, 11, 11, 22, 7915.6624, "B-", -4781.4051, 21.994437547, 0.141);
    static const Species Mg_22("Mg-22", "Mg", -2, 10, 12, 22, 7662.7645, "B-", -18601.0, 21.999570597, 0.17);
    static const Species Al_22("Al-22", "Al", -4, 9, 13, 22, 6782.0, "B-", -15439.0, 22.01954, 430.0);
    static const Species Si_22("Si-22", "Si", -6, 8, 14, 22, 6044.0, "B-", std::numeric_limits<double>::quiet_NaN(), 22.036114, 537.0);
    static const Species C_23("C-23", "C", 11, 17, 6, 23, 5077.0, "B-", 27450.0, 23.06889, 1070.0);
    static const Species N_23("N-23", "N", 9, 16, 7, 23, 6236.6721, "B-", 22099.0584, 23.039421, 451.5);
    static const Species O_23("O-23", "O", 7, 15, 8, 23, 7163.4856, "B-", 11336.1072, 23.015696686, 130.663);
    static const Species F_23("F-23", "F", 5, 14, 9, 23, 7622.3447, "B-", 8439.3084, 23.003526875, 35.77);
    static const Species Ne_23("Ne-23", "Ne", 3, 13, 10, 23, 7955.2561, "B-", 4375.8085, 22.994466905, 0.112);
    static const Species Na_23("Na-23", "Na", 1, 12, 11, 23, 8111.4936, "B-", -4056.179, 22.98976928195, 0.00194);
    static const Species Mg_23("Mg-23", "Mg", -1, 11, 12, 23, 7901.1229, "B-", -12221.7457, 22.994123768, 0.034);
    static const Species Al_23("Al-23", "Al", -3, 10, 13, 23, 7335.7275, "B-", -17202.0, 23.007244351, 0.37);
    static const Species Si_23("Si-23", "Si", -5, 9, 14, 23, 6554.0, "B-", std::numeric_limits<double>::quiet_NaN(), 23.025711, 537.0);
    static const Species N_24("N-24", "N", 10, 17, 7, 24, 5887.0, "B-", 28438.0, 24.05039, 430.0);
    static const Species O_24("O-24", "O", 8, 16, 8, 24, 7039.6855, "B-", 10955.8885, 24.019861, 177.0);
    static const Species F_24("F-24", "F", 6, 15, 9, 24, 7463.5831, "B-", 13496.1583, 24.00809937, 104.853);
    static const Species Ne_24("Ne-24", "Ne", 4, 14, 10, 24, 7993.3252, "B-", 2466.2583, 23.993610649, 0.55);
    static const Species Na_24("Na-24", "Na", 2, 13, 11, 24, 8063.4882, "B-", 5515.6774, 23.990963012, 0.017);
    static const Species Mg_24("Mg-24", "Mg", 0, 12, 12, 24, 8260.7103, "B-", -13884.766, 23.985041689, 0.013);
    static const Species Al_24("Al-24", "Al", -2, 11, 13, 24, 7649.5806, "B-", -10793.9978, 23.999947598, 0.244);
    static const Species Si_24("Si-24", "Si", -4, 10, 14, 24, 7167.2329, "B-", -23275.0, 24.01153543, 20.904);
    static const Species P_24("P-24", "P", -6, 9, 15, 24, 6165.0, "B-", std::numeric_limits<double>::quiet_NaN(), 24.036522, 537.0);
    static const Species N_25("N-25", "N", 11, 18, 7, 25, 5613.0, "B-", 28654.0, 25.0601, 540.0);
    static const Species O_25("O-25", "O", 9, 17, 8, 25, 6727.8058, "B-", 15994.8633, 25.029338919, 177.225);
    static const Species F_25("F-25", "F", 7, 16, 9, 25, 7336.3065, "B-", 13369.6698, 25.012167727, 103.535);
    static const Species Ne_25("Ne-25", "Ne", 5, 15, 10, 25, 7839.7994, "B-", 7322.3107, 24.997814797, 31.181);
    static const Species Na_25("Na-25", "Na", 3, 14, 11, 25, 8101.3979, "B-", 3834.9684, 24.989953974, 1.288);
    static const Species Mg_25("Mg-25", "Mg", 1, 13, 12, 25, 8223.5028, "B-", -4276.808, 24.985836966, 0.05);
    static const Species Al_25("Al-25", "Al", -1, 12, 13, 25, 8021.1366, "B-", -12743.2956, 24.990428308, 0.069);
    static const Species Si_25("Si-25", "Si", -3, 11, 14, 25, 7480.1109, "B-", -16363.0, 25.004108798, 10.735);
    static const Species P_25("P-25", "P", -5, 10, 15, 25, 6794.0, "B-", std::numeric_limits<double>::quiet_NaN(), 25.021675, 429.0);
    static const Species O_26("O-26", "O", 10, 18, 8, 26, 6497.479, "B-", 15986.3855, 26.037210155, 177.081);
    static const Species F_26("F-26", "F", 8, 17, 9, 26, 7082.2497, "B-", 18193.5414, 26.020048065, 114.898);
    static const Species Ne_26("Ne-26", "Ne", 6, 16, 10, 26, 7751.911, "B-", 7341.894, 26.000516496, 19.784);
    static const Species Na_26("Na-26", "Na", 4, 15, 11, 26, 8004.2013, "B-", 9353.7631, 25.992634649, 3.759);
    static const Species Mg_26("Mg-26", "Mg", 2, 14, 12, 26, 8333.8711, "B-", -4004.4042, 25.982592972, 0.031);
    static const Species Al_26("Al-26", "Al", 0, 13, 13, 26, 8149.7653, "B-", -5069.1361, 25.986891876, 0.071);
    static const Species Si_26("Si-26", "Si", -2, 12, 14, 26, 7924.7083, "B-", -18114.0, 25.992333818, 0.115);
    static const Species P_26("P-26", "P", -4, 11, 15, 26, 7198.0, "B-", -16707.0, 26.01178, 210.0);
    static const Species S_26("S-26", "S", -6, 10, 16, 26, 6525.0, "B-", std::numeric_limits<double>::quiet_NaN(), 26.029716, 644.0);
    static const Species O_27("O-27", "O", 11, 19, 8, 27, 6185.0, "B-", 19536.0, 27.047955, 537.0);
    static const Species F_27("F-27", "F", 9, 18, 9, 27, 6879.6662, "B-", 18082.568, 27.026981897, 129.037);
    static const Species Ne_27("Ne-27", "Ne", 7, 17, 10, 27, 7520.4151, "B-", 12568.7005, 27.007569462, 97.445);
    static const Species Na_27("Na-27", "Na", 5, 16, 11, 27, 7956.9467, "B-", 9068.8037, 26.994076408, 4.0);
    static const Species Mg_27("Mg-27", "Mg", 3, 15, 12, 27, 8263.8525, "B-", 2610.2694, 26.984340647, 0.05);
    static const Species Al_27("Al-27", "Al", 1, 14, 13, 27, 8331.5533, "B-", -4812.3583, 26.981538408, 0.05);
    static const Species Si_27("Si-27", "Si", -1, 13, 14, 27, 8124.342, "B-", -11725.473, 26.986704687, 0.115);
    static const Species P_27("P-27", "P", -3, 12, 15, 27, 7661.0894, "B-", -18150.0, 26.999292499, 9.662);
    static const Species S_27("S-27", "S", -5, 11, 16, 27, 6960.0, "B-", std::numeric_limits<double>::quiet_NaN(), 27.018777, 430.0);
    static const Species O_28("O-28", "O", 12, 20, 8, 28, 5988.0, "B-", 18676.0, 28.05591, 750.0);
    static const Species F_28("F-28", "F", 10, 19, 9, 28, 6626.8567, "B-", 22104.0579, 28.035860448, 129.198);
    static const Species Ne_28("Ne-28", "Ne", 8, 18, 10, 28, 7388.3463, "B-", 12288.0534, 28.012130767, 135.339);
    static const Species Na_28("Na-28", "Na", 6, 17, 11, 28, 7799.2644, "B-", 14031.6303, 27.998939, 11.0);
    static const Species Mg_28("Mg-28", "Mg", 4, 16, 12, 28, 8272.4531, "B-", 1830.774, 27.983875426, 0.28);
    static const Species Al_28("Al-28", "Al", 2, 15, 13, 28, 8309.8969, "B-", 4642.0776, 27.981910009, 0.052);
    static const Species Si_28("Si-28", "Si", 0, 14, 14, 28, 8447.7445, "B-", -14344.9407, 27.97692653442, 0.00055);
    static const Species P_28("P-28", "P", -2, 13, 15, 28, 7907.4842, "B-", -11221.0593, 27.99232646, 1.231);
    static const Species S_28("S-28", "S", -4, 12, 16, 28, 7478.7911, "B-", -24197.0, 28.004372762, 171.767);
    static const Species Cl_28("Cl-28", "Cl", -6, 11, 17, 28, 6587.0, "B-", std::numeric_limits<double>::quiet_NaN(), 28.030349, 537.0);
    static const Species F_29("F-29", "F", 11, 20, 9, 29, 6444.0314, "B-", 21750.3873, 29.043103, 564.0);
    static const Species Ne_29("Ne-29", "Ne", 9, 19, 10, 29, 7167.0673, "B-", 15719.8092, 29.019753, 160.5);
    static const Species Na_29("Na-29", "Na", 7, 18, 11, 29, 7682.1522, "B-", 13292.3538, 29.002877091, 7.876);
    static const Species Mg_29("Mg-29", "Mg", 5, 17, 12, 29, 8113.5317, "B-", 7595.4023, 28.988607163, 0.369);
    static const Species Al_29("Al-29", "Al", 3, 16, 13, 29, 8348.4647, "B-", 3687.3192, 28.980453164, 0.37);
    static const Species Si_29("Si-29", "Si", 1, 15, 14, 29, 8448.6361, "B-", -4942.2325, 28.97649466434, 0.0006);
    static const Species P_29("P-29", "P", -1, 14, 15, 29, 8251.2368, "B-", -13858.4257, 28.981800368000002, 0.385);
    static const Species S_29("S-29", "S", -3, 13, 16, 29, 7746.3826, "B-", -17117.0, 28.996678, 14.0);
    static const Species Cl_29("Cl-29", "Cl", -5, 12, 17, 29, 7129.0, "B-", -23947.0, 29.015053, 203.0);
    static const Species Ar_29("Ar-29", "Ar", -7, 11, 18, 29, 6276.0, "B-", std::numeric_limits<double>::quiet_NaN(), 29.040761, 471.0);
    static const Species F_30("F-30", "F", 12, 21, 9, 30, 6205.0, "B-", 25680.0, 30.052561, 537.0);
    static const Species Ne_30("Ne-30", "Ne", 10, 20, 10, 30, 7034.5317, "B-", 14805.4501, 30.024992235, 271.875);
    static const Species Na_30("Na-30", "Na", 8, 19, 11, 30, 7501.9685, "B-", 17356.0423, 30.009097931, 5.074);
    static const Species Mg_30("Mg-30", "Mg", 6, 18, 12, 30, 8054.425, "B-", 6982.744, 29.990465454, 1.39);
    static const Species Al_30("Al-30", "Al", 4, 17, 13, 30, 8261.1049, "B-", 8568.8459, 29.982969171, 2.077);
    static const Species Si_30("Si-30", "Si", 2, 16, 14, 30, 8520.6549, "B-", -4232.1065, 29.973770137, 0.023);
    static const Species P_30("P-30", "P", 0, 15, 15, 30, 8353.5064, "B-", -6141.6014, 29.97831349, 0.069);
    static const Species S_30("S-30", "S", -2, 14, 16, 30, 8122.7081, "B-", -18733.802, 29.98490677, 0.221);
    static const Species Cl_30("Cl-30", "Cl", -4, 13, 17, 30, 7472.1698, "B-", -17397.0, 30.005018333, 25.631);
    static const Species Ar_30("Ar-30", "Ar", -6, 12, 18, 30, 6866.0, "B-", std::numeric_limits<double>::quiet_NaN(), 30.023694, 192.0);
    static const Species F_31("F-31", "F", 13, 22, 9, 31, 6011.0, "B-", 25661.0, 31.061023, 574.0);
    static const Species Ne_31("Ne-31", "Ne", 11, 21, 10, 31, 6813.0902, "B-", 18935.5625, 31.033474816, 285.772);
    static const Species Na_31("Na-31", "Na", 9, 20, 11, 31, 7398.6778, "B-", 15368.1833, 31.013146654, 15.0);
    static const Species Mg_31("Mg-31", "Mg", 7, 19, 12, 31, 7869.1886, "B-", 11828.5569, 30.996648232, 3.3);
    static const Species Al_31("Al-31", "Al", 5, 18, 13, 31, 8225.518, "B-", 7998.3286, 30.983949754, 2.4);
    static const Species Si_31("Si-31", "Si", 3, 17, 14, 31, 8458.2916, "B-", 1491.5071, 30.975363196, 0.046);
    static const Species P_31("P-31", "P", 1, 16, 15, 31, 8481.1677, "B-", -5398.013, 30.97376199768, 0.0008);
    static const Species S_31("S-31", "S", -1, 15, 16, 31, 8281.8013, "B-", -12007.979, 30.979557002, 0.246);
    static const Species Cl_31("Cl-31", "Cl", -3, 14, 17, 31, 7869.2101, "B-", -18360.0, 30.992448097, 3.7);
    static const Species Ar_31("Ar-31", "Ar", -5, 13, 18, 31, 7252.0, "B-", -22935.0, 31.012158, 215.0);
    static const Species K_31("K-31", "K", -7, 12, 19, 31, 6487.0, "B-", std::numeric_limits<double>::quiet_NaN(), 31.03678, 322.0);
    static const Species Ne_32("Ne-32", "Ne", 12, 22, 10, 32, 6671.0, "B-", 18359.0, 32.03972, 540.0);
    static const Species Na_32("Na-32", "Na", 10, 21, 11, 32, 7219.8815, "B-", 19469.0523, 32.020011024, 40.0);
    static const Species Mg_32("Mg-32", "Mg", 8, 20, 12, 32, 7803.8411, "B-", 10270.4677, 31.999110138, 3.5);
    static const Species Al_32("Al-32", "Al", 6, 19, 13, 32, 8100.3449, "B-", 12978.3208, 31.988084338, 7.7);
    static const Species Si_32("Si-32", "Si", 4, 18, 14, 32, 8481.469, "B-", 227.1872, 31.974151538, 0.32);
    static const Species P_32("P-32", "P", 2, 17, 15, 32, 8464.1203, "B-", 1710.6608, 31.973907643, 0.042);
    static const Species S_32("S-32", "S", 0, 16, 16, 32, 8493.1301, "B-", -12680.8313, 31.97207117354, 0.00141);
    static const Species Cl_32("Cl-32", "Cl", -2, 15, 17, 32, 8072.4058, "B-", -11134.3536, 31.985684605, 0.603);
    static const Species Ar_32("Ar-32", "Ar", -4, 14, 18, 32, 7700.0089, "B-", -24190.0, 31.997637824, 1.9);
    static const Species K_32("K-32", "K", -6, 13, 19, 32, 6920.0, "B-", std::numeric_limits<double>::quiet_NaN(), 32.023607, 429.0);
    static const Species Ne_33("Ne-33", "Ne", 13, 23, 10, 33, 6436.0, "B-", 22350.0, 33.049523, 644.0);
    static const Species Na_33("Na-33", "Na", 11, 22, 11, 33, 7089.9262, "B-", 18817.24, 33.025529, 483.0);
    static const Species Mg_33("Mg-33", "Mg", 9, 21, 12, 33, 7636.4382, "B-", 13460.255, 33.005327862, 2.859);
    static const Species Al_33("Al-33", "Al", 7, 20, 13, 33, 8020.6172, "B-", 12016.946, 32.990877685, 7.5);
    static const Species Si_33("Si-33", "Si", 5, 19, 14, 33, 8361.0596, "B-", 5823.0223, 32.977976964, 0.75);
    static const Species P_33("P-33", "P", 3, 18, 15, 33, 8513.8073, "B-", 248.5079, 32.971725692, 1.17);
    static const Species S_33("S-33", "S", 1, 17, 16, 33, 8497.6304, "B-", -5582.5182, 32.97145890862, 0.00144);
    static const Species Cl_33("Cl-33", "Cl", -1, 16, 17, 33, 8304.7557, "B-", -11619.0452, 32.977451988, 0.419);
    static const Species Ar_33("Ar-33", "Ar", -3, 15, 18, 33, 7928.9559, "B-", -16925.0, 32.989925545, 0.43);
    static const Species K_33("K-33", "K", -5, 14, 19, 33, 7392.0, "B-", -23489.0, 33.008095, 215.0);
    static const Species Ca_33("Ca-33", "Ca", -7, 13, 20, 33, 6657.0, "B-", std::numeric_limits<double>::quiet_NaN(), 33.033312, 429.0);
    static const Species Ne_34("Ne-34", "Ne", 14, 24, 10, 34, 6287.0, "B-", 21161.0, 34.056728, 551.0);
    static const Species Na_34("Na-34", "Na", 12, 23, 11, 34, 6886.4377, "B-", 23356.7903, 34.03401, 643.5);
    static const Species Mg_34("Mg-34", "Mg", 10, 22, 12, 34, 7550.3919, "B-", 11320.9428, 34.008935455, 7.4);
    static const Species Al_34("Al-34", "Al", 8, 21, 13, 34, 7860.3506, "B-", 16994.0653, 33.996781924, 2.259);
    static const Species Si_34("Si-34", "Si", 6, 20, 14, 34, 8337.1659, "B-", 4557.0175, 33.978538045, 0.86);
    static const Species P_34("P-34", "P", 4, 19, 15, 34, 8448.1856, "B-", 5382.9879, 33.973645886, 0.87);
    static const Species S_34("S-34", "S", 2, 18, 16, 34, 8583.4986, "B-", -5491.6037, 33.967867011, 0.047);
    static const Species Cl_34("Cl-34", "Cl", 0, 17, 17, 34, 8398.9706, "B-", -6061.793, 33.97376249, 0.052);
    static const Species Ar_34("Ar-34", "Ar", -2, 16, 18, 34, 8197.6724, "B-", -17158.0, 33.980270092, 0.083);
    static const Species K_34("K-34", "K", -4, 15, 19, 34, 7670.0, "B-", -16110.0, 33.99869, 210.0);
    static const Species Ca_34("Ca-34", "Ca", -6, 14, 20, 34, 7173.0, "B-", std::numeric_limits<double>::quiet_NaN(), 34.015985, 322.0);
    static const Species Na_35("Na-35", "Na", 13, 24, 11, 35, 6745.0, "B-", 22192.0, 35.040614, 720.0);
    static const Species Mg_35("Mg-35", "Mg", 11, 23, 12, 35, 7356.2338, "B-", 15863.5156, 35.01679, 289.5);
    static const Species Al_35("Al-35", "Al", 9, 22, 13, 35, 7787.1243, "B-", 14167.7504, 34.999759816, 7.9);
    static const Species Si_35("Si-35", "Si", 7, 21, 14, 35, 8169.5644, "B-", 10466.3291, 34.984550111, 38.494);
    static const Species P_35("P-35", "P", 5, 20, 15, 35, 8446.2496, "B-", 3988.4006, 34.973314045, 2.003);
    static const Species S_35("S-35", "S", 3, 19, 16, 35, 8537.8511, "B-", 167.3218, 34.969032321, 0.043);
    static const Species Cl_35("Cl-35", "Cl", 1, 18, 17, 35, 8520.279, "B-", -5966.2429, 34.968852694, 0.038);
    static const Species Ar_35("Ar-35", "Ar", -1, 17, 18, 35, 8327.4621, "B-", -11874.3955, 34.975257719, 0.73);
    static const Species K_35("K-35", "K", -3, 16, 19, 35, 7965.8409, "B-", -16363.0, 34.988005406, 0.55);
    static const Species Ca_35("Ca-35", "Ca", -5, 15, 20, 35, 7476.0, "B-", -21910.0, 35.005572, 215.0);
    static const Species Sc_35("Sc-35", "Sc", -7, 14, 21, 35, 6828.0, "B-", std::numeric_limits<double>::quiet_NaN(), 35.029093, 429.0);
    static const Species Na_36("Na-36", "Na", 14, 25, 11, 36, 6557.0, "B-", 25523.0, 36.049279, 737.0);
    static const Species Mg_36("Mg-36", "Mg", 12, 24, 12, 36, 7244.4202, "B-", 14429.7751, 36.021879, 741.0);
    static const Species Al_36("Al-36", "Al", 10, 23, 13, 36, 7623.5154, "B-", 18386.5096, 36.006388, 160.5);
    static const Species Si_36("Si-36", "Si", 8, 22, 14, 36, 8112.5199, "B-", 7814.9194, 35.986649271, 77.077);
    static const Species P_36("P-36", "P", 6, 21, 15, 36, 8307.8692, "B-", 10413.0962, 35.97825961, 14.078);
    static const Species S_36("S-36", "S", 4, 20, 16, 36, 8575.39, "B-", -1142.1329, 35.967080692, 0.201);
    static const Species Cl_36("Cl-36", "Cl", 2, 19, 17, 36, 8521.9322, "B-", 709.5343, 35.968306822, 0.038);
    static const Species Ar_36("Ar-36", "Ar", 0, 18, 18, 36, 8519.9096, "B-", -12814.3607, 35.967545106, 0.028);
    static const Species K_36("K-36", "K", -2, 17, 19, 36, 8142.2233, "B-", -10966.0155, 35.981301887, 0.349);
    static const Species Ca_36("Ca-36", "Ca", -4, 16, 20, 36, 7815.8799, "B-", -22601.0, 35.993074388, 42.941);
    static const Species Sc_36("Sc-36", "Sc", -6, 15, 21, 36, 7166.0, "B-", std::numeric_limits<double>::quiet_NaN(), 36.017338, 322.0);
    static const Species Na_37("Na-37", "Na", 15, 26, 11, 37, 6403.0, "B-", 24923.0, 37.057042, 737.0);
    static const Species Mg_37("Mg-37", "Mg", 13, 25, 12, 37, 7055.1115, "B-", 18401.9132, 37.030286265, 750.35);
    static const Species Al_37("Al-37", "Al", 11, 24, 13, 37, 7531.316, "B-", 16381.0765, 37.010531, 193.5);
    static const Species Si_37("Si-37", "Si", 9, 23, 14, 37, 7952.9033, "B-", 12424.5003, 36.992945191, 122.179);
    static const Species P_37("P-37", "P", 7, 22, 15, 37, 8267.5561, "B-", 7900.4135, 36.979606942, 40.738);
    static const Species S_37("S-37", "S", 5, 21, 16, 37, 8459.9363, "B-", 4865.1258, 36.9711255, 0.212);
    static const Species Cl_37("Cl-37", "Cl", 3, 20, 17, 37, 8570.2816, "B-", -813.8729, 36.965902573, 0.055);
    static const Species Ar_37("Ar-37", "Ar", 1, 19, 18, 37, 8527.1406, "B-", -6147.4775, 36.966776301, 0.221);
    static const Species K_37("K-37", "K", -1, 18, 19, 37, 8339.848, "B-", -11664.1314, 36.97337589, 0.1);
    static const Species Ca_37("Ca-37", "Ca", -3, 17, 20, 37, 8003.4567, "B-", -16916.0, 36.985897849, 0.68);
    static const Species Sc_37("Sc-37", "Sc", -5, 16, 21, 37, 7525.0, "B-", -21390.0, 37.004058, 322.0);
    static const Species Ti_37("Ti-37", "Ti", -7, 15, 22, 37, 6926.0, "B-", std::numeric_limits<double>::quiet_NaN(), 37.027021, 429.0);
    static const Species Na_38("Na-38", "Na", 16, 27, 11, 38, 6216.0, "B-", 27831.0, 38.066458, 768.0);
    static const Species Mg_38("Mg-38", "Mg", 14, 26, 12, 38, 6928.0, "B-", 17604.0, 38.03658, 540.0);
    static const Species Al_38("Al-38", "Al", 12, 25, 13, 38, 7370.0, "B-", 20640.0, 38.017681, 161.0);
    static const Species Si_38("Si-38", "Si", 10, 24, 14, 38, 7892.8297, "B-", 10451.2656, 37.995523, 112.5);
    static const Species P_38("P-38", "P", 8, 23, 15, 38, 8147.2749, "B-", 12239.6512, 37.984303105, 77.918);
    static const Species S_38("S-38", "S", 6, 22, 16, 38, 8448.7829, "B-", 2936.9, 37.9711633, 7.699);
    static const Species Cl_38("Cl-38", "Cl", 4, 21, 17, 38, 8505.4817, "B-", 4916.7109, 37.968010408, 0.105);
    static const Species Ar_38("Ar-38", "Ar", 2, 20, 18, 38, 8614.2807, "B-", -5914.0671, 37.962732102, 0.209);
    static const Species K_38("K-38", "K", 0, 19, 19, 38, 8438.0593, "B-", -6742.2563, 37.969081114, 0.209);
    static const Species Ca_38("Ca-38", "Ca", -2, 18, 20, 38, 8240.0434, "B-", -17809.0, 37.976319223, 0.208);
    static const Species Sc_38("Sc-38", "Sc", -4, 17, 21, 38, 7751.0, "B-", -15619.0, 37.995438, 215.0);
    static const Species Ti_38("Ti-38", "Ti", -6, 16, 22, 38, 7319.0, "B-", std::numeric_limits<double>::quiet_NaN(), 38.012206, 322.0);
    static const Species Na_39("Na-39", "Na", 17, 28, 11, 39, 6056.0, "B-", 27201.0, 39.075123, 797.0);
    static const Species Mg_39("Mg-39", "Mg", 15, 27, 12, 39, 6734.0, "B-", 21286.0, 39.045921, 551.0);
    static const Species Al_39("Al-39", "Al", 13, 26, 13, 39, 7260.0, "B-", 19169.0, 39.02307, 322.0);
    static const Species Si_39("Si-39", "Si", 11, 25, 14, 39, 7730.9793, "B-", 15094.9874, 39.002491, 145.5);
    static const Species P_39("P-39", "P", 9, 24, 15, 39, 8097.9701, "B-", 10388.0361, 38.986285865, 120.929);
    static const Species S_39("S-39", "S", 7, 23, 16, 39, 8344.2698, "B-", 6637.5463, 38.97513385, 53.677);
    static const Species Cl_39("Cl-39", "Cl", 5, 22, 17, 39, 8494.4032, "B-", 3441.9774, 38.968008151, 1.859);
    static const Species Ar_39("Ar-39", "Ar", 3, 21, 18, 39, 8562.5988, "B-", 565.0, 38.964313037, 5.367);
    static const Species K_39("K-39", "K", 1, 20, 19, 39, 8557.0258, "B-", -6524.4888, 38.96370648482, 0.00489);
    static const Species Ca_39("Ca-39", "Ca", -1, 19, 20, 39, 8369.6711, "B-", -13109.9804, 38.970710811000004, 0.64);
    static const Species Sc_39("Sc-39", "Sc", -3, 18, 21, 39, 8013.4575, "B-", -16673.0, 38.984784953, 25.765);
    static const Species Ti_39("Ti-39", "Ti", -5, 17, 22, 39, 7566.0, "B-", -20070.0, 39.002684, 215.0);
    static const Species V_39("V-39", "V", -7, 16, 23, 39, 7031.0, "B-", std::numeric_limits<double>::quiet_NaN(), 39.02423, 429.0);
    static const Species Mg_40("Mg-40", "Mg", 16, 28, 12, 40, 6598.0, "B-", 20729.0, 40.053194, 537.0);
    static const Species Al_40("Al-40", "Al", 14, 27, 13, 40, 7097.0, "B-", 23154.0, 40.03094, 322.0);
    static const Species Si_40("Si-40", "Si", 12, 26, 14, 40, 7655.8247, "B-", 13806.0653, 40.006083641, 130.962);
    static const Species P_40("P-40", "P", 10, 25, 15, 40, 7981.4177, "B-", 14698.6603, 39.991262221, 89.755);
    static const Species S_40("S-40", "S", 8, 24, 16, 40, 8329.3255, "B-", 4719.9687, 39.975482561, 4.274);
    static const Species Cl_40("Cl-40", "Cl", 6, 23, 17, 40, 8427.766, "B-", 7482.0816, 39.970415466, 34.423);
    static const Species Ar_40("Ar-40", "Ar", 4, 22, 18, 40, 8595.2594, "B-", -1504.4031, 39.96238312204, 0.00234);
    static const Species K_40("K-40", "K", 2, 21, 19, 40, 8538.0907, "B-", 1310.9051, 39.963998165, 0.06);
    static const Species Ca_40("Ca-40", "Ca", 0, 20, 20, 40, 8551.3046, "B-", -14323.0493, 39.96259085, 0.022);
    static const Species Sc_40("Sc-40", "Sc", -2, 19, 21, 40, 8173.6697, "B-", -11529.9139, 39.977967275, 3.036);
    static const Species Ti_40("Ti-40", "Ti", -4, 18, 22, 40, 7865.8632, "B-", -21463.0, 39.990345146, 73.262);
    static const Species V_40("V-40", "V", -6, 17, 23, 40, 7310.0, "B-", std::numeric_limits<double>::quiet_NaN(), 40.013387, 322.0);
    static const Species Mg_41("Mg-41", "Mg", 17, 29, 12, 41, 6425.0, "B-", 23510.0, 41.062373, 537.0);
    static const Species Al_41("Al-41", "Al", 15, 28, 13, 41, 6980.0, "B-", 21390.0, 41.037134, 429.0);
    static const Species Si_41("Si-41", "Si", 13, 27, 14, 41, 7482.0, "B-", 18180.0, 41.014171, 322.0);
    static const Species P_41("P-41", "P", 11, 26, 15, 41, 7906.5513, "B-", 14028.8125, 40.994654, 129.0);
    static const Species S_41("S-41", "S", 9, 25, 16, 41, 8229.6358, "B-", 8298.6116, 40.979593451, 4.4);
    static const Species Cl_41("Cl-41", "Cl", 7, 24, 17, 41, 8412.9593, "B-", 5760.318, 40.970684525, 73.777);
    static const Species Ar_41("Ar-41", "Ar", 5, 23, 18, 41, 8534.3733, "B-", 2492.0392, 40.96450057, 0.372);
    static const Species K_41("K-41", "K", 3, 22, 19, 41, 8576.0731, "B-", -421.6406, 40.96182525611, 0.00403);
    static const Species Ca_41("Ca-41", "Ca", 1, 21, 20, 41, 8546.7075, "B-", -6495.5482, 40.962277905, 0.147);
    static const Species Sc_41("Sc-41", "Sc", -1, 20, 21, 41, 8369.1979, "B-", -12944.8214, 40.969251163, 0.083);
    static const Species Ti_41("Ti-41", "Ti", -3, 19, 22, 41, 8034.3889, "B-", -16008.0, 40.983148, 30.0);
    static const Species V_41("V-41", "V", -5, 18, 23, 41, 7625.0, "B-", -20100.0, 41.000333, 215.0);
    static const Species Cr_41("Cr-41", "Cr", -7, 17, 24, 41, 7116.0, "B-", std::numeric_limits<double>::quiet_NaN(), 41.021911, 429.0);
    static const Species Al_42("Al-42", "Al", 16, 29, 13, 42, 6829.0, "B-", 25150.0, 42.045078, 537.0);
    static const Species Si_42("Si-42", "Si", 14, 28, 14, 42, 7410.0, "B-", 15748.0, 42.018078, 322.0);
    static const Species P_42("P-42", "P", 12, 27, 15, 42, 7765.9122, "B-", 18729.5899, 42.00117214, 101.996);
    static const Species S_42("S-42", "S", 10, 26, 16, 42, 8193.2275, "B-", 7194.0221, 41.9810651, 3.0);
    static const Species Cl_42("Cl-42", "Cl", 8, 25, 17, 42, 8345.8864, "B-", 9590.9082, 41.973342, 64.0);
    static const Species Ar_42("Ar-42", "Ar", 6, 24, 18, 42, 8555.6141, "B-", 599.3527, 41.963045737, 6.2);
    static const Species K_42("K-42", "K", 4, 23, 19, 42, 8551.2571, "B-", 3525.2626, 41.962402305, 0.113);
    static const Species Ca_42("Ca-42", "Ca", 2, 22, 20, 42, 8616.5646, "B-", -6426.2904, 41.95861778, 0.159);
    static const Species Sc_42("Sc-42", "Sc", 0, 21, 21, 42, 8444.9303, "B-", -7016.6496, 41.965516686, 0.165);
    static const Species Ti_42("Ti-42", "Ti", -2, 20, 22, 42, 8259.24, "B-", -17485.0, 41.973049369, 0.289);
    static const Species V_42("V-42", "V", -4, 19, 23, 42, 7824.0, "B-", -14679.0, 41.99182, 210.0);
    static const Species Cr_42("Cr-42", "Cr", -6, 18, 24, 42, 7456.0, "B-", std::numeric_limits<double>::quiet_NaN(), 42.007579, 322.0);
    static const Species Al_43("Al-43", "Al", 17, 30, 13, 43, 6712.0, "B-", 23940.0, 43.05182, 644.0);
    static const Species Si_43("Si-43", "Si", 15, 29, 14, 43, 7251.0, "B-", 19289.0, 43.026119, 429.0);
    static const Species P_43("P-43", "P", 13, 28, 15, 43, 7681.0, "B-", 17236.0, 43.005411, 322.0);
    static const Species S_43("S-43", "S", 11, 27, 16, 43, 8063.8276, "B-", 11964.05, 42.986907635, 5.335);
    static const Species Cl_43("Cl-43", "Cl", 9, 26, 17, 43, 8323.8672, "B-", 7850.3002, 42.9740637, 66.407);
    static const Species Ar_43("Ar-43", "Ar", 7, 25, 18, 43, 8488.2382, "B-", 4565.5836, 42.965636056, 5.7);
    static const Species K_43("K-43", "K", 5, 24, 19, 43, 8576.2204, "B-", 1833.4783, 42.960734701, 0.44);
    static const Species Ca_43("Ca-43", "Ca", 3, 23, 20, 43, 8600.6653, "B-", -2220.7227, 42.958766381, 0.244);
    static const Species Sc_43("Sc-43", "Sc", 1, 22, 21, 43, 8530.8265, "B-", -6872.5591, 42.961150425, 1.999);
    static const Species Ti_43("Ti-43", "Ti", -1, 21, 22, 43, 8352.8054, "B-", -11399.2333, 42.96852842, 6.139);
    static const Species V_43("V-43", "V", -3, 20, 23, 43, 8069.5129, "B-", -15946.0, 42.980766, 46.0);
    static const Species Cr_43("Cr-43", "Cr", -5, 19, 24, 43, 7680.0, "B-", -19340.0, 42.997885, 215.0);
    static const Species Mn_43("Mn-43", "Mn", -7, 18, 25, 43, 7213.0, "B-", std::numeric_limits<double>::quiet_NaN(), 43.018647, 429.0);
    static const Species Si_44("Si-44", "Si", 16, 30, 14, 44, 7156.0, "B-", 18200.0, 44.031466, 537.0);
    static const Species P_44("P-44", "P", 14, 29, 15, 44, 7552.0, "B-", 20314.0, 44.011927, 429.0);
    static const Species S_44("S-44", "S", 12, 28, 16, 44, 7996.0154, "B-", 11274.7373, 43.990118846, 5.6);
    static const Species Cl_44("Cl-44", "Cl", 10, 27, 17, 44, 8234.4788, "B-", 12194.2869, 43.978014918, 91.859);
    static const Species Ar_44("Ar-44", "Ar", 8, 26, 18, 44, 8493.8411, "B-", 3108.2375, 43.964923814, 1.7);
    static const Species K_44("K-44", "K", 6, 25, 19, 44, 8546.7023, "B-", 5687.2319, 43.961586984, 0.45);
    static const Species Ca_44("Ca-44", "Ca", 4, 24, 20, 44, 8658.1769, "B-", -3652.6948, 43.955481489, 0.348);
    static const Species Sc_44("Sc-44", "Sc", 2, 23, 21, 44, 8557.3805, "B-", -267.4488, 43.959402818, 1.884);
    static const Species Ti_44("Ti-44", "Ti", 0, 22, 22, 44, 8533.5215, "B-", -13740.5071, 43.959689936, 0.751);
    static const Species V_44("V-44", "V", -2, 21, 23, 44, 8203.4567, "B-", -10386.1805, 43.974440977, 7.799);
    static const Species Cr_44("Cr-44", "Cr", -4, 20, 24, 44, 7949.6265, "B-", -20882.0, 43.985591, 55.0);
    static const Species Mn_44("Mn-44", "Mn", -6, 19, 25, 44, 7457.0, "B-", std::numeric_limits<double>::quiet_NaN(), 44.008009, 322.0);
    static const Species Si_45("Si-45", "Si", 17, 31, 14, 45, 7004.0, "B-", 21130.0, 45.039818, 644.0);
    static const Species P_45("P-45", "P", 15, 30, 15, 45, 7456.0, "B-", 19301.0, 45.017134, 537.0);
    static const Species S_45("S-45", "S", 13, 29, 16, 45, 7867.0, "B-", 14922.0, 44.996414, 322.0);
    static const Species Cl_45("Cl-45", "Cl", 11, 28, 17, 45, 8181.5991, "B-", 11508.2568, 44.980394353, 146.177);
    static const Species Ar_45("Ar-45", "Ar", 9, 27, 18, 45, 8419.9526, "B-", 6844.8422, 44.968039731, 0.55);
    static const Species K_45("K-45", "K", 7, 26, 19, 45, 8554.6747, "B-", 4196.5868, 44.960691491, 0.56);
    static const Species Ca_45("Ca-45", "Ca", 5, 25, 20, 45, 8630.5467, "B-", 260.091, 44.95618627, 0.392);
    static const Species Sc_45("Sc-45", "Sc", 3, 24, 21, 45, 8618.941, "B-", -2062.0551, 44.955907051, 0.712);
    static const Species Ti_45("Ti-45", "Ti", 1, 23, 22, 45, 8555.7321, "B-", -7123.8247, 44.958120758, 0.897);
    static const Species V_45("V-45", "V", -1, 22, 23, 45, 8380.0394, "B-", -12371.64, 44.965768498, 0.926);
    static const Species Cr_45("Cr-45", "Cr", -3, 21, 24, 45, 8087.7286, "B-", -14535.0, 44.97905, 38.0);
    static const Species Mn_45("Mn-45", "Mn", -5, 20, 25, 45, 7747.0, "B-", -19388.0, 44.994654, 322.0);
    static const Species Fe_45("Fe-45", "Fe", -7, 19, 26, 45, 7299.0, "B-", std::numeric_limits<double>::quiet_NaN(), 45.015467, 304.0);
    static const Species P_46("P-46", "P", 16, 31, 15, 46, 7320.0, "B-", 22200.0, 46.02452, 537.0);
    static const Species S_46("S-46", "S", 14, 30, 16, 46, 7785.0, "B-", 14375.0, 46.000687, 429.0);
    static const Species Cl_46("Cl-46", "Cl", 12, 29, 17, 46, 8080.7757, "B-", 16036.3062, 45.985254926, 104.4);
    static const Species Ar_46("Ar-46", "Ar", 10, 28, 18, 46, 8412.3835, "B-", 5642.6746, 45.968039244, 2.5);
    static const Species K_46("K-46", "K", 8, 27, 19, 46, 8518.0428, "B-", 7725.6802, 45.961981584, 0.78);
    static const Species Ca_46("Ca-46", "Ca", 6, 26, 20, 46, 8668.9848, "B-", -1377.9665, 45.953687726, 2.398);
    static const Species Sc_46("Sc-46", "Sc", 4, 25, 21, 46, 8622.0215, "B-", 2366.626, 45.955167034, 0.72);
    static const Species Ti_46("Ti-46", "Ti", 2, 24, 22, 46, 8656.4623, "B-", -7052.3723, 45.952626356, 0.097);
    static const Species V_46("V-46", "V", 0, 23, 23, 46, 8486.1423, "B-", -7604.3264, 45.960197389, 0.143);
    static const Species Cr_46("Cr-46", "Cr", -2, 22, 24, 46, 8303.8233, "B-", -17053.8226, 45.968360969, 12.295);
    static const Species Mn_46("Mn-46", "Mn", -4, 21, 25, 46, 7916.0805, "B-", -13628.0, 45.986669, 93.0);
    static const Species Fe_46("Fe-46", "Fe", -6, 20, 26, 46, 7603.0, "B-", std::numeric_limits<double>::quiet_NaN(), 46.001299, 322.0);
    static const Species P_47("P-47", "P", 17, 32, 15, 47, 7209.0, "B-", 21610.0, 47.030929, 644.0);
    static const Species S_47("S-47", "S", 15, 31, 16, 47, 7652.0, "B-", 16781.0, 47.00773, 429.0);
    static const Species Cl_47("Cl-47", "Cl", 13, 30, 17, 47, 7992.0, "B-", 15787.0, 46.989715, 215.0);
    static const Species Ar_47("Ar-47", "Ar", 11, 29, 18, 47, 8311.425, "B-", 10344.7078, 46.972767112, 1.3);
    static const Species K_47("K-47", "K", 9, 28, 19, 47, 8514.8795, "B-", 6632.6837, 46.961661612, 1.5);
    static const Species Ca_47("Ca-47", "Ca", 7, 27, 20, 47, 8639.3548, "B-", 1992.177, 46.954541134, 2.384);
    static const Species Sc_47("Sc-47", "Sc", 5, 26, 21, 47, 8665.0958, "B-", 600.7694, 46.952402444, 2.072);
    static const Species Ti_47("Ti-47", "Ti", 3, 25, 22, 47, 8661.2325, "B-", -2930.5422, 46.951757491, 0.085);
    static const Species V_47("V-47", "V", 1, 24, 23, 47, 8582.2348, "B-", -7443.9769, 46.954903558, 0.118);
    static const Species Cr_47("Cr-47", "Cr", -1, 23, 24, 47, 8407.2067, "B-", -11996.7167, 46.962894995, 5.578);
    static const Species Mn_47("Mn-47", "Mn", -3, 22, 25, 47, 8135.3117, "B-", -15437.0, 46.975774, 34.0);
    static const Species Fe_47("Fe-47", "Fe", -5, 21, 26, 47, 7790.0, "B-", -17750.0, 46.992346, 537.0);
    static const Species Co_47("Co-47", "Co", -7, 20, 27, 47, 7396.0, "B-", std::numeric_limits<double>::quiet_NaN(), 47.011401, 644.0);
    static const Species S_48("S-48", "S", 16, 32, 16, 48, 7552.0, "B-", 16670.0, 48.013301, 537.0);
    static const Species Cl_48("Cl-48", "Cl", 14, 31, 17, 48, 7883.0, "B-", 18075.0, 47.995405, 537.0);
    static const Species Ar_48("Ar-48", "Ar", 12, 30, 18, 48, 8243.6656, "B-", 9929.555, 47.976001, 18.0);
    static const Species K_48("K-48", "K", 10, 29, 19, 48, 8434.2324, "B-", 11940.3857, 47.965341184, 0.83);
    static const Species Ca_48("Ca-48", "Ca", 8, 28, 20, 48, 8666.6916, "B-", 279.2155, 47.952522654, 0.018);
    static const Species Sc_48("Sc-48", "Sc", 6, 27, 21, 48, 8656.2097, "B-", 3988.8685, 47.952222903, 5.313);
    static const Species Ti_48("Ti-48", "Ti", 4, 26, 22, 48, 8723.0122, "B-", -4014.9467, 47.947940677, 0.079);
    static const Species V_48("V-48", "V", 2, 25, 23, 48, 8623.0686, "B-", -1656.6918, 47.9522509, 1.043);
    static const Species Cr_48("Cr-48", "Cr", 0, 24, 24, 48, 8572.2553, "B-", -13524.6692, 47.954029431, 7.848);
    static const Species Mn_48("Mn-48", "Mn", -2, 23, 25, 48, 8274.1924, "B-", -11288.0685, 47.96854876, 7.191);
    static const Species Fe_48("Fe-48", "Fe", -4, 22, 26, 48, 8022.7254, "B-", -19738.0, 47.980667, 99.0);
    static const Species Co_48("Co-48", "Co", -6, 21, 27, 48, 7595.0, "B-", -16448.0, 48.001857, 537.0);
    static const Species Ni_48("Ni-48", "Ni", -8, 20, 28, 48, 7236.0, "B-", std::numeric_limits<double>::quiet_NaN(), 48.019515, 455.0);
    static const Species S_49("S-49", "S", 17, 33, 16, 49, 7400.0, "B-", 19652.0, 49.021891, 626.0);
    static const Species Cl_49("Cl-49", "Cl", 15, 32, 17, 49, 7785.0, "B-", 17800.0, 49.000794, 429.0);
    static const Species Ar_49("Ar-49", "Ar", 13, 31, 18, 49, 8132.0, "B-", 12551.0, 48.981685, 429.0);
    static const Species K_49("K-49", "K", 11, 30, 19, 49, 8372.2753, "B-", 11688.5069, 48.968210753, 0.86);
    static const Species Ca_49("Ca-49", "Ca", 9, 29, 20, 49, 8594.85, "B-", 5262.4445, 48.955662625, 0.19);
    static const Species Sc_49("Sc-49", "Sc", 7, 28, 21, 49, 8686.2805, "B-", 2001.5652, 48.950013159, 2.434);
    static const Species Ti_49("Ti-49", "Ti", 5, 27, 22, 49, 8711.1625, "B-", -601.8555, 48.947864391, 0.084);
    static const Species V_49("V-49", "V", 3, 26, 23, 49, 8682.9135, "B-", -2629.8047, 48.948510509, 0.884);
    static const Species Cr_49("Cr-49", "Cr", 1, 25, 24, 49, 8613.2777, "B-", -7712.4265, 48.95133372, 2.363);
    static const Species Mn_49("Mn-49", "Mn", -1, 24, 25, 49, 8439.915, "B-", -12869.1957, 48.95961335, 2.377);
    static const Species Fe_49("Fe-49", "Fe", -3, 23, 26, 49, 8161.3121, "B-", -14971.0, 48.973429, 26.0);
    static const Species Co_49("Co-49", "Co", -5, 22, 27, 49, 7840.0, "B-", -18309.0, 48.989501, 537.0);
    static const Species Ni_49("Ni-49", "Ni", -7, 21, 28, 49, 7450.0, "B-", std::numeric_limits<double>::quiet_NaN(), 49.009157, 644.0);
    static const Species Cl_50("Cl-50", "Cl", 16, 33, 17, 50, 7651.0, "B-", 20930.0, 50.008266, 429.0);
    static const Species Ar_50("Ar-50", "Ar", 14, 32, 18, 50, 8054.0, "B-", 12498.0, 49.985797, 537.0);
    static const Species K_50("K-50", "K", 12, 31, 19, 50, 8288.5833, "B-", 13861.3774, 49.972380015, 8.3);
    static const Species Ca_50("Ca-50", "Ca", 10, 30, 20, 50, 8550.1639, "B-", 4947.8903, 49.957499215, 1.7);
    static const Species Sc_50("Sc-50", "Sc", 8, 29, 21, 50, 8633.4747, "B-", 6894.747, 49.952187437, 2.7);
    static const Species Ti_50("Ti-50", "Ti", 6, 28, 22, 50, 8755.7227, "B-", -2208.6274, 49.944785622, 0.088);
    static const Species V_50("V-50", "V", 4, 27, 23, 50, 8695.9032, "B-", 1038.124, 49.947156681, 0.099);
    static const Species Cr_50("Cr-50", "Cr", 2, 26, 24, 50, 8701.0188, "B-", -7634.4776, 49.946042209, 0.1);
    static const Species Mn_50("Mn-50", "Mn", 0, 25, 25, 50, 8532.6823, "B-", -8150.4267, 49.954238157, 0.123);
    static const Species Fe_50("Fe-50", "Fe", -2, 24, 26, 50, 8354.0268, "B-", -16887.0566, 49.962988, 9.0);
    static const Species Co_50("Co-50", "Co", -4, 23, 27, 50, 8000.6387, "B-", -14130.0, 49.981117, 135.0);
    static const Species Ni_50("Ni-50", "Ni", -6, 22, 28, 50, 7702.0, "B-", std::numeric_limits<double>::quiet_NaN(), 49.996286, 537.0);
    static const Species Cl_51("Cl-51", "Cl", 17, 34, 17, 51, 7530.0, "B-", 20780.0, 51.015341, 751.0);
    static const Species Ar_51("Ar-51", "Ar", 15, 33, 18, 51, 7922.0, "B-", 16026.0, 50.993033, 429.0);
    static const Species K_51("K-51", "K", 13, 32, 19, 51, 8221.335, "B-", 13816.8529, 50.975828664, 14.0);
    static const Species Ca_51("Ca-51", "Ca", 11, 31, 20, 51, 8476.9135, "B-", 6918.0432, 50.960995663, 0.56);
    static const Species Sc_51("Sc-51", "Sc", 9, 30, 21, 51, 8597.2213, "B-", 6482.6122, 50.953568838, 2.7);
    static const Species Ti_51("Ti-51", "Ti", 7, 29, 22, 51, 8708.9912, "B-", 2470.1402, 50.946609468, 0.519);
    static const Species V_51("V-51", "V", 5, 28, 23, 51, 8742.0852, "B-", -752.3907, 50.943957664, 0.104);
    static const Species Cr_51("Cr-51", "Cr", 3, 27, 24, 51, 8711.9923, "B-", -3207.4893, 50.944765388, 0.178);
    static const Species Mn_51("Mn-51", "Mn", 1, 26, 25, 51, 8633.7602, "B-", -8054.04, 50.94820877, 0.326);
    static const Species Fe_51("Fe-51", "Fe", -1, 25, 26, 51, 8460.4977, "B-", -12847.0389, 50.956855137, 1.501);
    static const Species Co_51("Co-51", "Co", -3, 24, 27, 51, 8193.2549, "B-", -15692.0, 50.970647, 52.0);
    static const Species Ni_51("Ni-51", "Ni", -5, 23, 28, 51, 7870.0, "B-", std::numeric_limits<double>::quiet_NaN(), 50.987493, 537.0);
    static const Species Cl_52("Cl-52", "Cl", 18, 35, 17, 52, 7386.0, "B-", 23739.0, 52.024004, 751.0);
    static const Species Ar_52("Ar-52", "Ar", 16, 34, 18, 52, 7827.0, "B-", 15758.0, 51.998519, 644.0);
    static const Species K_52("K-52", "K", 14, 33, 19, 52, 8115.0303, "B-", 17128.6431, 51.981602, 36.0);
    static const Species Ca_52("Ca-52", "Ca", 12, 32, 20, 52, 8429.3821, "B-", 6257.2889, 51.963213646, 0.72);
    static const Species Sc_52("Sc-52", "Sc", 10, 31, 21, 52, 8534.6695, "B-", 8954.1372, 51.95649617, 3.3);
    static const Species Ti_52("Ti-52", "Ti", 8, 30, 22, 52, 8691.8193, "B-", 1965.334, 51.946883509, 2.948);
    static const Species V_52("V-52", "V", 6, 29, 23, 52, 8714.569, "B-", 3976.4763, 51.944773636, 0.17);
    static const Species Cr_52("Cr-52", "Cr", 4, 28, 24, 52, 8775.9946, "B-", -4708.1214, 51.940504714, 0.12);
    static const Species Mn_52("Mn-52", "Mn", 2, 27, 25, 52, 8670.4087, "B-", -2379.2912, 51.94555909, 0.138);
    static const Species Fe_52("Fe-52", "Fe", 0, 26, 26, 52, 8609.6079, "B-", -13988.1167, 51.948113364, 0.192);
    static const Species Co_52("Co-52", "Co", -2, 25, 27, 52, 8325.5606, "B-", -11784.123, 51.963130224, 5.669);
    static const Species Ni_52("Ni-52", "Ni", -4, 24, 28, 52, 8083.8977, "B-", -20680.0, 51.975781, 89.0);
    static const Species Cu_52("Cu-52", "Cu", -6, 23, 29, 52, 7671.0, "B-", std::numeric_limits<double>::quiet_NaN(), 51.997982, 644.0);
    static const Species Ar_53("Ar-53", "Ar", 17, 35, 18, 53, 7677.0, "B-", 19086.0, 53.00729, 750.0);
    static const Species K_53("K-53", "K", 15, 34, 19, 53, 8022.8488, "B-", 17091.9853, 52.9868, 120.0);
    static const Species Ca_53("Ca-53", "Ca", 13, 33, 20, 53, 8330.5778, "B-", 9381.8471, 52.968451, 47.0);
    static const Species Sc_53("Sc-53", "Sc", 11, 32, 21, 53, 8492.8325, "B-", 8111.8785, 52.958379173, 19.0);
    static const Species Ti_53("Ti-53", "Ti", 9, 31, 22, 53, 8631.1256, "B-", 4970.2419, 52.949670714, 3.1);
    static const Species V_53("V-53", "V", 7, 30, 23, 53, 8710.1425, "B-", 3435.9426, 52.94433494, 3.331);
    static const Species Cr_53("Cr-53", "Cr", 5, 29, 24, 53, 8760.2103, "B-", -597.2679, 52.940646304, 0.124);
    static const Species Mn_53("Mn-53", "Mn", 3, 28, 25, 53, 8734.1798, "B-", -3742.8664, 52.941287497, 0.371);
    static const Species Fe_53("Fe-53", "Fe", 1, 27, 26, 53, 8648.7985, "B-", -8288.1073, 52.945305629, 1.792);
    static const Species Co_53("Co-53", "Co", -1, 26, 27, 53, 8477.6578, "B-", -13028.5485, 52.954203278, 1.854);
    static const Species Ni_53("Ni-53", "Ni", -3, 25, 28, 53, 8217.0749, "B-", -16491.0, 52.96819, 27.0);
    static const Species Cu_53("Cu-53", "Cu", -5, 24, 29, 53, 7891.0, "B-", std::numeric_limits<double>::quiet_NaN(), 52.985894, 537.0);
    static const Species Ar_54("Ar-54", "Ar", 18, 36, 18, 54, 7578.0, "B-", 17710.0, 54.013484, 859.0);
    static const Species K_54("K-54", "K", 16, 35, 19, 54, 7891.0, "B-", 20010.0, 53.994471, 429.0);
    static const Species Ca_54("Ca-54", "Ca", 14, 34, 20, 54, 8247.4967, "B-", 9277.3463, 53.972989, 52.0);
    static const Species Sc_54("Sc-54", "Sc", 12, 33, 21, 54, 8404.8115, "B-", 11305.8789, 53.963029359000004, 15.0);
    static const Species Ti_54("Ti-54", "Ti", 10, 32, 22, 54, 8599.6917, "B-", 4154.4548, 53.950892, 17.0);
    static const Species V_54("V-54", "V", 8, 31, 23, 54, 8662.1382, "B-", 7037.1118, 53.946432009, 12.001);
    static const Species Cr_54("Cr-54", "Cr", 6, 30, 24, 54, 8777.9672, "B-", -1377.1325, 53.938877359, 0.142);
    static const Species Mn_54("Mn-54", "Mn", 4, 29, 25, 54, 8737.9768, "B-", 696.3688, 53.940355772, 1.08);
    static const Species Fe_54("Fe-54", "Fe", 2, 28, 26, 54, 8736.3846, "B-", -8244.5478, 53.939608189, 0.368);
    static const Species Co_54("Co-54", "Co", 0, 27, 27, 54, 8569.2199, "B-", -8731.7558, 53.948459075, 0.38);
    static const Species Ni_54("Ni-54", "Ni", -2, 26, 28, 54, 8393.0328, "B-", -18038.0, 53.957833, 5.0);
    static const Species Cu_54("Cu-54", "Cu", -4, 25, 29, 54, 8045.0, "B-", -15538.0, 53.977198, 429.0);
    static const Species Zn_54("Zn-54", "Zn", -6, 24, 30, 54, 7742.0, "B-", std::numeric_limits<double>::quiet_NaN(), 53.993879, 232.0);
    static const Species K_55("K-55", "K", 17, 36, 19, 55, 7792.0, "B-", 19121.0, 55.000505, 537.0);
    static const Species Ca_55("Ca-55", "Ca", 15, 35, 20, 55, 8125.926, "B-", 12191.7329, 54.979978, 172.0);
    static const Species Sc_55("Sc-55", "Sc", 13, 34, 21, 55, 8333.3694, "B-", 10990.3608, 54.966889637, 67.0);
    static const Species Ti_55("Ti-55", "Ti", 11, 33, 22, 55, 8518.9696, "B-", 7292.6673, 54.955091, 31.0);
    static const Species V_55("V-55", "V", 9, 32, 23, 55, 8637.3391, "B-", 5985.1877, 54.947262, 29.0);
    static const Species Cr_55("Cr-55", "Cr", 7, 31, 24, 55, 8731.9362, "B-", 2602.2183, 54.940836637, 0.245);
    static const Species Mn_55("Mn-55", "Mn", 5, 30, 25, 55, 8765.0247, "B-", -231.1204, 54.93804304, 0.279);
    static const Species Fe_55("Fe-55", "Fe", 3, 29, 26, 55, 8746.5981, "B-", -3451.4254, 54.938291158, 0.33);
    static const Species Co_55("Co-55", "Co", 1, 28, 27, 55, 8669.6204, "B-", -8694.035, 54.941996416, 0.434);
    static const Species Ni_55("Ni-55", "Ni", -1, 27, 28, 55, 8497.3225, "B-", -13700.5585, 54.951329846, 0.757);
    static const Species Cu_55("Cu-55", "Cu", -3, 26, 29, 55, 8233.997, "B-", -17366.0, 54.966038, 167.0);
    static const Species Zn_55("Zn-55", "Zn", -5, 25, 30, 55, 7904.0, "B-", std::numeric_limits<double>::quiet_NaN(), 54.984681, 429.0);
    static const Species K_56("K-56", "K", 18, 37, 19, 56, 7663.0, "B-", 21491.0, 56.008567, 644.0);
    static const Species Ca_56("Ca-56", "Ca", 16, 36, 20, 56, 8033.1654, "B-", 12005.458, 55.985496, 268.0);
    static const Species Sc_56("Sc-56", "Sc", 14, 35, 21, 56, 8233.5781, "B-", 13907.147, 55.972607611, 278.761);
    static const Species Ti_56("Ti-56", "Ti", 12, 34, 22, 56, 8467.9495, "B-", 6760.4051, 55.957677675, 107.569);
    static const Species V_56("V-56", "V", 10, 33, 23, 56, 8574.7006, "B-", 9101.727, 55.950420082, 188.819);
    static const Species Cr_56("Cr-56", "Cr", 8, 32, 24, 56, 8723.2609, "B-", 1626.5384, 55.940648977, 0.62);
    static const Species Mn_56("Mn-56", "Mn", 6, 31, 25, 56, 8738.3358, "B-", 3695.4973, 55.938902816, 0.314);
    static const Species Fe_56("Fe-56", "Fe", 4, 30, 26, 56, 8790.3563, "B-", -4566.6455, 55.934935537, 0.287);
    static const Species Co_56("Co-56", "Co", 2, 29, 27, 56, 8694.8386, "B-", -2132.8689, 55.939838032, 0.51);
    static const Species Ni_56("Ni-56", "Ni", 0, 28, 28, 56, 8642.7811, "B-", -15277.9163, 55.942127761, 0.428);
    static const Species Cu_56("Cu-56", "Cu", -2, 27, 29, 56, 8355.9907, "B-", -13240.0, 55.958529278, 6.864);
    static const Species Zn_56("Zn-56", "Zn", -4, 26, 30, 56, 8106.0, "B-", -21550.0, 55.972743, 429.0);
    static const Species Ga_56("Ga-56", "Ga", -6, 25, 31, 56, 7707.0, "B-", std::numeric_limits<double>::quiet_NaN(), 55.995878, 537.0);
    static const Species K_57("K-57", "K", 19, 38, 19, 57, 7563.0, "B-", 20689.0, 57.015169, 644.0);
    static const Species Ca_57("Ca-57", "Ca", 17, 37, 20, 57, 7912.0, "B-", 14820.0, 56.992958, 429.0);
    static const Species Sc_57("Sc-57", "Sc", 15, 36, 21, 57, 8158.1666, "B-", 13022.1957, 56.977048, 193.0);
    static const Species Ti_57("Ti-57", "Ti", 13, 35, 22, 57, 8372.9008, "B-", 10033.2148, 56.963068098, 221.02);
    static const Species V_57("V-57", "V", 11, 34, 23, 57, 8535.1967, "B-", 8089.9214, 56.952297, 91.0);
    static const Species Cr_57("Cr-57", "Cr", 9, 33, 24, 57, 8663.3998, "B-", 4961.2946, 56.943612112, 2.0);
    static const Species Mn_57("Mn-57", "Mn", 7, 32, 25, 57, 8736.7146, "B-", 2695.7375, 56.938285944, 1.615);
    static const Species Fe_57("Fe-57", "Fe", 5, 31, 26, 57, 8770.2829, "B-", -836.3589, 56.93539195, 0.287);
    static const Species Co_57("Co-57", "Co", 3, 30, 27, 57, 8741.8845, "B-", -3261.697, 56.936289819, 0.553);
    static const Species Ni_57("Ni-57", "Ni", 1, 29, 28, 57, 8670.9364, "B-", -8774.9466, 56.939791394, 0.608);
    static const Species Cu_57("Cu-57", "Cu", -1, 28, 29, 57, 8503.2646, "B-", -14759.0, 56.949211686, 0.537);
    static const Species Zn_57("Zn-57", "Zn", -3, 27, 30, 57, 8231.0, "B-", -17140.0, 56.965056, 215.0);
    static const Species Ga_57("Ga-57", "Ga", -5, 26, 31, 57, 7916.0, "B-", std::numeric_limits<double>::quiet_NaN(), 56.983457, 429.0);
    static const Species K_58("K-58", "K", 20, 39, 19, 58, 7437.0, "B-", 23461.0, 58.023543, 751.0);
    static const Species Ca_58("Ca-58", "Ca", 18, 38, 20, 58, 7828.0, "B-", 13949.0, 57.998357, 537.0);
    static const Species Sc_58("Sc-58", "Sc", 16, 37, 21, 58, 8054.9436, "B-", 15438.0995, 57.983382, 204.0);
    static const Species Ti_58("Ti-58", "Ti", 14, 36, 22, 58, 8307.629, "B-", 9512.9152, 57.966808519, 196.823);
    static const Species V_58("V-58", "V", 12, 35, 23, 58, 8458.156, "B-", 11561.224, 57.956595985, 102.862);
    static const Species Cr_58("Cr-58", "Cr", 10, 34, 24, 58, 8643.9987, "B-", 3835.7607, 57.944184501, 3.2);
    static const Species Mn_58("Mn-58", "Mn", 8, 33, 25, 58, 8696.6438, "B-", 6327.7027, 57.940066643, 2.9);
    static const Species Fe_58("Fe-58", "Fe", 6, 32, 26, 58, 8792.2534, "B-", -2307.9785, 57.933273575, 0.339);
    static const Species Co_58("Co-58", "Co", 4, 31, 27, 58, 8738.9719, "B-", 381.5789, 57.935751292, 1.237);
    static const Species Ni_58("Ni-58", "Ni", 2, 30, 28, 58, 8732.0621, "B-", -8561.0204, 57.93534165, 0.374);
    static const Species Cu_58("Cu-58", "Cu", 0, 29, 29, 58, 8570.9696, "B-", -9368.9796, 57.944532283, 0.604);
    static const Species Zn_58("Zn-58", "Zn", -2, 28, 30, 58, 8395.9467, "B-", -18759.0, 57.954590296, 53.678);
    static const Species Ga_58("Ga-58", "Ga", -4, 27, 31, 58, 8059.0, "B-", -15960.0, 57.974728999999996, 322.0);
    static const Species Ge_58("Ge-58", "Ge", -6, 26, 32, 58, 7770.0, "B-", std::numeric_limits<double>::quiet_NaN(), 57.991863, 537.0);
    static const Species K_59("K-59", "K", 21, 40, 19, 59, 7332.0, "B-", 22940.0, 59.030864, 859.0);
    static const Species Ca_59("Ca-59", "Ca", 19, 39, 20, 59, 7708.0, "B-", 16639.0, 59.006237, 644.0);
    static const Species Sc_59("Sc-59", "Sc", 17, 38, 21, 59, 7976.4073, "B-", 15050.0, 58.988374, 268.0);
    static const Species Ti_59("Ti-59", "Ti", 15, 37, 22, 59, 8218.0, "B-", 11731.0, 58.972217, 322.0);
    static const Species V_59("V-59", "V", 13, 36, 23, 59, 8403.8034, "B-", 10505.3137, 58.959623343, 147.505);
    static const Species Cr_59("Cr-59", "Cr", 11, 35, 24, 59, 8568.5995, "B-", 7409.3977, 58.948345426, 0.72);
    static const Species Mn_59("Mn-59", "Mn", 9, 34, 25, 59, 8680.9224, "B-", 5139.629, 58.940391111, 2.5);
    static const Species Fe_59("Fe-59", "Fe", 7, 33, 26, 59, 8754.7746, "B-", 1564.8804, 58.934873492, 0.354);
    static const Species Co_59("Co-59", "Co", 5, 32, 27, 59, 8768.0379, "B-", -1073.005, 58.933193524, 0.426);
    static const Species Ni_59("Ni-59", "Ni", 3, 31, 28, 59, 8736.5912, "B-", -4798.3786, 58.934345442, 0.376);
    static const Species Cu_59("Cu-59", "Cu", 1, 30, 29, 59, 8642.0027, "B-", -9142.776, 58.939496713, 0.566);
    static const Species Zn_59("Zn-59", "Zn", -1, 29, 30, 59, 8473.7802, "B-", -13456.0, 58.949311886, 0.814);
    static const Species Ga_59("Ga-59", "Ga", -3, 28, 31, 59, 8232.0, "B-", -17390.0, 58.963757, 183.0);
    static const Species Ge_59("Ge-59", "Ge", -5, 27, 32, 59, 7924.0, "B-", std::numeric_limits<double>::quiet_NaN(), 58.982426, 429.0);
    static const Species Ca_60("Ca-60", "Ca", 20, 40, 20, 60, 7627.0, "B-", 15550.0, 60.011809, 751.0);
    static const Species Sc_60("Sc-60", "Sc", 18, 39, 21, 60, 7873.0, "B-", 17549.0, 59.995115, 537.0);
    static const Species Ti_60("Ti-60", "Ti", 16, 38, 22, 60, 8152.7858, "B-", 10987.704, 59.976275, 258.0);
    static const Species V_60("V-60", "V", 14, 37, 23, 60, 8322.8751, "B-", 13821.0986, 59.964479215, 195.327);
    static const Species Cr_60("Cr-60", "Cr", 12, 36, 24, 60, 8540.1876, "B-", 6059.4457, 59.949641656, 1.2);
    static const Species Mn_60("Mn-60", "Mn", 10, 35, 25, 60, 8628.1392, "B-", 8445.2283, 59.943136574, 2.5);
    static const Species Fe_60("Fe-60", "Fe", 8, 34, 26, 60, 8755.8539, "B-", 237.2633, 59.934070249, 3.656);
    static const Species Co_60("Co-60", "Co", 6, 33, 27, 60, 8746.7692, "B-", 2822.8058, 59.933815536, 0.433);
    static const Species Ni_60("Ni-60", "Ni", 4, 32, 28, 60, 8780.7769, "B-", -6127.981, 59.930785129, 0.378);
    static const Species Cu_60("Cu-60", "Cu", 2, 31, 29, 60, 8665.6047, "B-", -4170.7922, 59.937363787, 1.731);
    static const Species Zn_60("Zn-60", "Zn", 0, 30, 30, 60, 8583.0524, "B-", -14584.0, 59.941841317, 0.588);
    static const Species Ga_60("Ga-60", "Ga", -2, 29, 31, 60, 8327.0, "B-", -12060.0, 59.957498, 215.0);
    static const Species Ge_60("Ge-60", "Ge", -4, 28, 32, 60, 8113.0, "B-", -21890.0, 59.970445, 322.0);
    static const Species As_60("As-60", "As", -6, 27, 33, 60, 7735.0, "B-", std::numeric_limits<double>::quiet_NaN(), 59.993945, 429.0);
    static const Species Ca_61("Ca-61", "Ca", 21, 41, 20, 61, 7503.0, "B-", 18510.0, 61.020408, 859.0);
    static const Species Sc_61("Sc-61", "Sc", 19, 40, 21, 61, 7794.0, "B-", 16870.0, 61.000537, 644.0);
    static const Species Ti_61("Ti-61", "Ti", 17, 39, 22, 61, 8058.0, "B-", 13807.0, 60.982426, 322.0);
    static const Species V_61("V-61", "V", 15, 38, 23, 61, 8271.0417, "B-", 12319.3807, 60.967603529, 252.196);
    static const Species Cr_61("Cr-61", "Cr", 13, 37, 24, 61, 8460.1734, "B-", 9245.6277, 60.95437813, 2.0);
    static const Species Mn_61("Mn-61", "Mn", 11, 36, 25, 61, 8598.9157, "B-", 7178.373, 60.944452541, 2.5);
    static const Species Fe_61("Fe-61", "Fe", 9, 35, 26, 61, 8703.7686, "B-", 3977.6759, 60.936746241, 2.8);
    static const Species Co_61("Co-61", "Co", 7, 34, 27, 61, 8756.151, "B-", 1323.8504, 60.932476031, 0.901);
    static const Species Ni_61("Ni-61", "Ni", 5, 33, 28, 61, 8765.0281, "B-", -2237.9663, 60.931054819, 0.381);
    static const Species Cu_61("Cu-61", "Cu", 3, 32, 29, 61, 8715.5148, "B-", -5635.1565, 60.933457375, 1.02);
    static const Species Zn_61("Zn-61", "Zn", 1, 31, 30, 61, 8610.3098, "B-", -9214.2438, 60.939506964, 17.068);
    static const Species Ga_61("Ga-61", "Ga", -1, 30, 31, 61, 8446.4313, "B-", -13345.0, 60.949398861, 40.787);
    static const Species Ge_61("Ge-61", "Ge", -3, 29, 32, 61, 8215.0, "B-", -16590.0, 60.963725, 322.0);
    static const Species As_61("As-61", "As", -5, 28, 33, 61, 7930.0, "B-", std::numeric_limits<double>::quiet_NaN(), 60.981535, 322.0);
    static const Species Sc_62("Sc-62", "Sc", 20, 41, 21, 62, 7688.0, "B-", 19510.0, 62.007848, 644.0);
    static const Species Ti_62("Ti-62", "Ti", 18, 40, 22, 62, 7990.0, "B-", 13013.0, 61.986903, 429.0);
    static const Species V_62("V-62", "V", 16, 39, 23, 62, 8187.7565, "B-", 15639.4478, 61.972932556, 283.723);
    static const Species Cr_62("Cr-62", "Cr", 14, 38, 24, 62, 8427.3871, "B-", 7671.3532, 61.95614292, 3.7);
    static const Species Mn_62("Mn-62", "Mn", 12, 37, 25, 62, 8538.5002, "B-", 10354.092, 61.947907384, 7.023);
    static const Species Fe_62("Fe-62", "Fe", 10, 36, 26, 62, 8692.8831, "B-", 2546.3427, 61.936791809, 3.0);
    static const Species Co_62("Co-62", "Co", 8, 35, 27, 62, 8721.3347, "B-", 5322.0404, 61.934058198, 19.94);
    static const Species Ni_62("Ni-62", "Ni", 6, 34, 28, 62, 8794.5555, "B-", -3958.8965, 61.928344753, 0.455);
    static const Species Cu_62("Cu-62", "Cu", 4, 33, 29, 62, 8718.0839, "B-", -1619.4548, 61.932594803, 0.683);
    static const Species Zn_62("Zn-62", "Zn", 2, 32, 30, 62, 8679.3451, "B-", -9181.0666, 61.934333359, 0.66);
    static const Species Ga_62("Ga-62", "Ga", 0, 31, 31, 62, 8518.6449, "B-", -9847.0, 61.944189639, 0.684);
    static const Species Ge_62("Ge-62", "Ge", -2, 30, 32, 62, 8347.0, "B-", -17720.0, 61.954761, 150.0);
    static const Species As_62("As-62", "As", -4, 29, 33, 62, 8049.0, "B-", std::numeric_limits<double>::quiet_NaN(), 61.973784, 322.0);
    static const Species Sc_63("Sc-63", "Sc", 21, 42, 21, 63, 7603.0, "B-", 18930.0, 63.014031, 751.0);
    static const Species Ti_63("Ti-63", "Ti", 19, 41, 22, 63, 7891.0, "B-", 15880.0, 62.993709, 537.0);
    static const Species V_63("V-63", "V", 17, 40, 23, 63, 8130.7809, "B-", 14438.1586, 62.976661, 365.0);
    static const Species Cr_63("Cr-63", "Cr", 15, 39, 24, 63, 8347.5398, "B-", 10708.7611, 62.961161, 78.0);
    static const Species Mn_63("Mn-63", "Mn", 13, 38, 25, 63, 8505.102, "B-", 8748.5685, 62.949664672, 4.0);
    static const Species Fe_63("Fe-63", "Fe", 11, 37, 26, 63, 8631.5499, "B-", 6215.9238, 62.940272698, 4.618);
    static const Species Co_63("Co-63", "Co", 9, 36, 27, 63, 8717.7972, "B-", 3661.3385, 62.93359963, 19.941);
    static const Species Ni_63("Ni-63", "Ni", 7, 35, 28, 63, 8763.4955, "B-", 66.9768, 62.929669021, 0.457);
    static const Species Cu_63("Cu-63", "Cu", 5, 34, 29, 63, 8752.1404, "B-", -3366.4392, 62.929597119, 0.457);
    static const Species Zn_63("Zn-63", "Zn", 3, 33, 30, 63, 8686.2866, "B-", -5666.3294, 62.93321114, 1.674);
    static const Species Ga_63("Ga-63", "Ga", 1, 32, 31, 63, 8583.9267, "B-", -9625.8787, 62.939294194, 1.4);
    static const Species Ge_63("Ge-63", "Ge", -1, 31, 32, 63, 8418.7167, "B-", -13421.0, 62.949628, 40.0);
    static const Species As_63("As-63", "As", -3, 30, 33, 63, 8193.0, "B-", -16650.0, 62.964036, 215.0);
    static const Species Se_63("Se-63", "Se", -5, 29, 34, 63, 7917.0, "B-", std::numeric_limits<double>::quiet_NaN(), 62.981911, 537.0);
    static const Species Ti_64("Ti-64", "Ti", 20, 42, 22, 64, 7826.0, "B-", 14840.0, 63.998411, 644.0);
    static const Species V_64("V-64", "V", 18, 41, 23, 64, 8045.0, "B-", 17320.0, 63.98248, 429.0);
    static const Species Cr_64("Cr-64", "Cr", 16, 40, 24, 64, 8303.5626, "B-", 9349.0622, 63.963886, 322.0);
    static const Species Mn_64("Mn-64", "Mn", 14, 39, 25, 64, 8437.4175, "B-", 11980.5117, 63.953849369, 3.8);
    static const Species Fe_64("Fe-64", "Fe", 12, 38, 26, 64, 8612.3888, "B-", 4822.8898, 63.940987761, 5.386);
    static const Species Co_64("Co-64", "Co", 10, 37, 27, 64, 8675.5223, "B-", 7306.5921, 63.935810176, 21.476);
    static const Species Ni_64("Ni-64", "Ni", 8, 36, 28, 64, 8777.4637, "B-", -1674.6156, 63.927966228, 0.497);
    static const Species Cu_64("Cu-64", "Cu", 6, 35, 29, 64, 8739.0736, "B-", 579.5996, 63.929764001, 0.458);
    static const Species Zn_64("Zn-64", "Zn", 4, 34, 30, 64, 8735.9057, "B-", -7171.1912, 63.929141776, 0.69);
    static const Species Ga_64("Ga-64", "Ga", 2, 33, 31, 64, 8611.6317, "B-", -4517.3237, 63.936840366, 1.533);
    static const Species Ge_64("Ge-64", "Ge", 0, 32, 32, 64, 8528.8243, "B-", -14783.0, 63.941689912, 4.0);
    static const Species As_64("As-64", "As", -2, 31, 33, 64, 8286.0, "B-", -12673.0, 63.95756, 218.0);
    static const Species Se_64("Se-64", "Se", -4, 30, 34, 64, 8075.0, "B-", std::numeric_limits<double>::quiet_NaN(), 63.971165, 537.0);
    static const Species Ti_65("Ti-65", "Ti", 21, 43, 22, 65, 7726.0, "B-", 17320.0, 65.005593, 751.0);
    static const Species V_65("V-65", "V", 19, 42, 23, 65, 7981.0, "B-", 16200.0, 64.986999, 537.0);
    static const Species Cr_65("Cr-65", "Cr", 17, 41, 24, 65, 8218.0, "B-", 12657.0, 64.969608, 215.0);
    static const Species Mn_65("Mn-65", "Mn", 15, 40, 25, 65, 8400.6822, "B-", 10250.5576, 64.956019749, 4.0);
    static const Species Fe_65("Fe-65", "Fe", 13, 39, 26, 65, 8546.347, "B-", 7967.3036, 64.945015323, 5.487);
    static const Species Co_65("Co-65", "Co", 11, 38, 27, 65, 8656.8848, "B-", 5940.5911, 64.936462071, 2.235);
    static const Species Ni_65("Ni-65", "Ni", 9, 37, 28, 65, 8736.2424, "B-", 2137.8808, 64.930084585, 0.518);
    static const Species Cu_65("Cu-65", "Cu", 7, 36, 29, 65, 8757.0967, "B-", -1351.6527, 64.927789476, 0.69);
    static const Species Zn_65("Zn-65", "Zn", 5, 35, 30, 65, 8724.266, "B-", -3254.538, 64.929240534, 0.693);
    static const Species Ga_65("Ga-65", "Ga", 3, 34, 31, 65, 8662.1601, "B-", -6179.2631, 64.932734424, 0.849);
    static const Species Ge_65("Ge-65", "Ge", 1, 33, 32, 65, 8555.0584, "B-", -9541.167, 64.939368136, 2.323);
    static const Species As_65("As-65", "As", -1, 32, 33, 65, 8396.2351, "B-", -13917.0, 64.949611, 91.0);
    static const Species Se_65("Se-65", "Se", -3, 31, 34, 65, 8170.0, "B-", -16529.0, 64.964552, 322.0);
    static const Species Br_65("Br-65", "Br", -5, 30, 35, 65, 7904.0, "B-", std::numeric_limits<double>::quiet_NaN(), 64.982297, 537.0);
    static const Species V_66("V-66", "V", 20, 43, 23, 66, 7894.0, "B-", 18840.0, 65.993237, 537.0);
    static const Species Cr_66("Cr-66", "Cr", 18, 42, 24, 66, 8168.0, "B-", 11610.0, 65.973011, 322.0);
    static const Species Mn_66("Mn-66", "Mn", 16, 41, 25, 66, 8331.7986, "B-", 13317.4543, 65.960546833, 12.0);
    static const Species Fe_66("Fe-66", "Fe", 14, 40, 26, 66, 8521.7245, "B-", 6340.6944, 65.946249958, 4.4);
    static const Species Co_66("Co-66", "Co", 12, 39, 27, 66, 8605.9419, "B-", 9597.7522, 65.939442943, 15.0);
    static const Species Ni_66("Ni-66", "Ni", 10, 38, 28, 66, 8739.5086, "B-", 251.9958, 65.929139333, 1.5);
    static const Species Cu_66("Cu-66", "Cu", 8, 37, 29, 66, 8731.473, "B-", 2640.9396, 65.928868804, 0.696);
    static const Species Zn_66("Zn-66", "Zn", 6, 36, 30, 66, 8759.6335, "B-", -5175.5, 65.926033639, 0.798);
    static const Species Ga_66("Ga-66", "Ga", 4, 35, 31, 66, 8669.3631, "B-", -2116.6879, 65.931589766, 1.172);
    static const Species Ge_66("Ge-66", "Ge", 2, 34, 32, 66, 8625.4383, "B-", -9581.957, 65.933862124, 2.577);
    static const Species As_66("As-66", "As", 0, 33, 33, 66, 8468.4034, "B-", -10365.0, 65.944148778, 6.1);
    static const Species Se_66("Se-66", "Se", -2, 32, 34, 66, 8300.0, "B-", -18091.0, 65.955276, 215.0);
    static const Species Br_66("Br-66", "Br", -4, 31, 35, 66, 8014.0, "B-", std::numeric_limits<double>::quiet_NaN(), 65.974697, 429.0);
    static const Species V_67("V-67", "V", 21, 44, 23, 67, 7829.0, "B-", 17526.0, 66.998128, 644.0);
    static const Species Cr_67("Cr-67", "Cr", 19, 43, 24, 67, 8079.0, "B-", 14311.0, 66.979313, 429.0);
    static const Species Mn_67("Mn-67", "Mn", 17, 42, 25, 67, 8281.0, "B-", 12128.0, 66.96395, 215.0);
    static const Species Fe_67("Fe-67", "Fe", 15, 41, 26, 67, 8449.9359, "B-", 9613.3678, 66.95093, 4.1);
    static const Species Co_67("Co-67", "Co", 13, 40, 27, 67, 8581.7422, "B-", 8420.9047, 66.940609625, 6.917);
    static const Species Ni_67("Ni-67", "Ni", 11, 39, 28, 67, 8695.7505, "B-", 3576.8654, 66.931569413, 3.1);
    static const Species Cu_67("Cu-67", "Cu", 9, 38, 29, 67, 8737.4597, "B-", 560.8226, 66.92772949, 0.957);
    static const Species Zn_67("Zn-67", "Zn", 7, 37, 30, 67, 8734.1534, "B-", -1001.2201, 66.927127422, 0.81);
    static const Species Ga_67("Ga-67", "Ga", 5, 36, 31, 67, 8707.533, "B-", -4205.438, 66.928202276, 1.262);
    static const Species Ge_67("Ge-67", "Ge", 3, 35, 32, 67, 8633.0884, "B-", -6086.4858, 66.932716999, 4.636);
    static const Species As_67("As-67", "As", 1, 34, 33, 67, 8530.5685, "B-", -10006.9381, 66.93925111, 0.475);
    static const Species Se_67("Se-67", "Se", -1, 33, 34, 67, 8369.5344, "B-", -14051.0, 66.949994, 72.0);
    static const Species Br_67("Br-67", "Br", -3, 32, 35, 67, 8148.0, "B-", -16978.0, 66.965078, 322.0);
    static const Species Kr_67("Kr-67", "Kr", -5, 31, 36, 67, 7883.0, "B-", std::numeric_limits<double>::quiet_NaN(), 66.983305, 455.0);
    static const Species Cr_68("Cr-68", "Cr", 20, 44, 24, 68, 8026.0, "B-", 13230.0, 67.983156, 537.0);
    static const Species Mn_68("Mn-68", "Mn", 18, 43, 25, 68, 8209.0, "B-", 14977.0, 67.968953, 322.0);
    static const Species Fe_68("Fe-68", "Fe", 16, 42, 26, 68, 8418.0, "B-", 7746.0, 67.952875, 207.0);
    static const Species Co_68("Co-68", "Co", 14, 41, 27, 68, 8520.1301, "B-", 11821.2318, 67.944559401, 4.142);
    static const Species Ni_68("Ni-68", "Ni", 12, 40, 28, 68, 8682.4667, "B-", 2103.2205, 67.931868787, 3.2);
    static const Species Cu_68("Cu-68", "Cu", 10, 39, 29, 68, 8701.8913, "B-", 4440.1115, 67.929610887, 1.7);
    static const Species Zn_68("Zn-68", "Zn", 8, 38, 30, 68, 8755.682, "B-", -2921.1, 67.924844232, 0.835);
    static const Species Ga_68("Ga-68", "Ga", 6, 37, 31, 68, 8701.2195, "B-", -107.2555, 67.927980161, 1.535);
    static const Species Ge_68("Ge-68", "Ge", 4, 36, 32, 68, 8688.1371, "B-", -8084.2715, 67.928095305, 2.014);
    static const Species As_68("As-68", "As", 2, 35, 33, 68, 8557.7457, "B-", -4705.0786, 67.936774127, 1.981);
    static const Species Se_68("Se-68", "Se", 0, 34, 34, 68, 8477.0482, "B-", -15398.0, 67.941825236, 0.532);
    static const Species Br_68("Br-68", "Br", -2, 33, 35, 68, 8239.0, "B-", -13165.0, 67.958356, 278.0);
    static const Species Kr_68("Kr-68", "Kr", -4, 32, 36, 68, 8034.0, "B-", std::numeric_limits<double>::quiet_NaN(), 67.972489, 537.0);
    static const Species Cr_69("Cr-69", "Cr", 21, 45, 24, 69, 7939.0, "B-", 15730.0, 68.989662, 537.0);
    static const Species Mn_69("Mn-69", "Mn", 19, 44, 25, 69, 8155.0, "B-", 13839.0, 68.972775, 429.0);
    static const Species Fe_69("Fe-69", "Fe", 17, 43, 26, 69, 8345.0, "B-", 11186.0, 68.957918, 215.0);
    static const Species Co_69("Co-69", "Co", 15, 42, 27, 69, 8495.4062, "B-", 9593.2084, 68.945909, 92.0);
    static const Species Ni_69("Ni-69", "Ni", 13, 41, 28, 69, 8623.0998, "B-", 5757.565, 68.935610267, 4.0);
    static const Species Cu_69("Cu-69", "Cu", 11, 40, 29, 69, 8695.2044, "B-", 2681.6854, 68.929429267, 1.5);
    static const Species Zn_69("Zn-69", "Zn", 9, 39, 30, 69, 8722.7311, "B-", 909.9134, 68.92655036, 0.853);
    static const Species Ga_69("Ga-69", "Ga", 7, 38, 31, 69, 8724.5798, "B-", -2227.1455, 68.925573528, 1.285);
    static const Species Ge_69("Ge-69", "Ge", 5, 37, 32, 69, 8680.964, "B-", -3988.4927, 68.927964467, 1.414);
    static const Species As_69("As-69", "As", 3, 36, 33, 69, 8611.8214, "B-", -6677.4672, 68.932246289, 34.352);
    static const Species Se_69("Se-69", "Se", 1, 35, 34, 69, 8503.7082, "B-", -10175.2364, 68.939414845, 1.599);
    static const Species Br_69("Br-69", "Br", -1, 34, 35, 69, 8344.9026, "B-", -14119.0, 68.95033841, 45.091);
    static const Species Kr_69("Kr-69", "Kr", -3, 33, 36, 69, 8129.0, "B-", std::numeric_limits<double>::quiet_NaN(), 68.965496, 322.0);
    static const Species Cr_70("Cr-70", "Cr", 22, 46, 24, 70, 7884.0, "B-", 14810.0, 69.993945, 644.0);
    static const Species Mn_70("Mn-70", "Mn", 20, 45, 25, 70, 8084.0, "B-", 16440.0, 69.978046, 537.0);
    static const Species Fe_70("Fe-70", "Fe", 18, 44, 26, 70, 8308.0, "B-", 9635.0, 69.960397, 322.0);
    static const Species Co_70("Co-70", "Co", 16, 43, 27, 70, 8434.198, "B-", 12688.9049, 69.9500534, 11.8);
    static const Species Ni_70("Ni-70", "Ni", 14, 42, 28, 70, 8604.2917, "B-", 3762.5123, 69.9364313, 2.301);
    static const Species Cu_70("Cu-70", "Cu", 12, 41, 29, 70, 8646.8655, "B-", 6588.3675, 69.932392078, 1.161);
    static const Species Zn_70("Zn-70", "Zn", 10, 40, 30, 70, 8729.8086, "B-", -654.5979, 69.925319175, 2.058);
    static const Species Ga_70("Ga-70", "Ga", 8, 39, 31, 70, 8709.2808, "B-", 1651.8861, 69.926021914, 1.289);
    static const Species Ge_70("Ge-70", "Ge", 6, 38, 32, 70, 8721.7028, "B-", -6228.063, 69.924248542, 0.88);
    static const Species As_70("As-70", "As", 4, 37, 33, 70, 8621.5541, "B-", -2404.0737, 69.930934642, 1.5);
    static const Species Se_70("Se-70", "Se", 2, 36, 34, 70, 8576.0338, "B-", -10504.2727, 69.933515521, 1.7);
    static const Species Br_70("Br-70", "Br", 0, 35, 35, 70, 8414.7964, "B-", -10325.0, 69.944792321, 16.0);
    static const Species Kr_70("Kr-70", "Kr", -2, 34, 36, 70, 8256.0, "B-", std::numeric_limits<double>::quiet_NaN(), 69.955877, 215.0);
    static const Species Mn_71("Mn-71", "Mn", 21, 46, 25, 71, 8030.0, "B-", 15310.0, 70.982158, 537.0);
    static const Species Fe_71("Fe-71", "Fe", 19, 45, 26, 71, 8235.0, "B-", 12440.0, 70.965722, 429.0);
    static const Species Co_71("Co-71", "Co", 17, 44, 27, 71, 8398.7344, "B-", 11036.3053, 70.952366923, 499.23);
    static const Species Ni_71("Ni-71", "Ni", 15, 43, 28, 71, 8543.1564, "B-", 7304.8989, 70.940518962, 2.401);
    static const Species Cu_71("Cu-71", "Cu", 13, 42, 29, 71, 8635.0233, "B-", 4617.6517, 70.932676831, 1.6);
    static const Species Zn_71("Zn-71", "Zn", 11, 41, 30, 71, 8689.0417, "B-", 2810.3405, 70.927719578, 2.849);
    static const Species Ga_71("Ga-71", "Ga", 9, 40, 31, 71, 8717.605, "B-", -232.4698, 70.924702554, 0.87);
    static const Species Ge_71("Ge-71", "Ge", 7, 39, 32, 71, 8703.3118, "B-", -2013.4, 70.92495212, 0.874);
    static const Species As_71("As-71", "As", 5, 38, 33, 71, 8663.935, "B-", -4746.742, 70.927113594, 4.469);
    static const Species Se_71("Se-71", "Se", 3, 37, 34, 71, 8586.0606, "B-", -6644.0883, 70.932209431, 3.0);
    static const Species Br_71("Br-71", "Br", 1, 36, 35, 71, 8481.4629, "B-", -10175.2155, 70.939342153, 5.799);
    static const Species Kr_71("Kr-71", "Kr", -1, 35, 36, 71, 8327.131, "B-", -14037.0, 70.950265695, 138.238);
    static const Species Rb_71("Rb-71", "Rb", -3, 34, 37, 71, 8118.0, "B-", std::numeric_limits<double>::quiet_NaN(), 70.965335, 429.0);
    static const Species Mn_72("Mn-72", "Mn", 22, 47, 25, 72, 7955.0, "B-", 18080.0, 71.988009, 644.0);
    static const Species Fe_72("Fe-72", "Fe", 20, 46, 26, 72, 8195.0, "B-", 11050.0, 71.968599, 537.0);
    static const Species Co_72("Co-72", "Co", 18, 45, 27, 72, 8338.0, "B-", 13926.0, 71.956736, 322.0);
    static const Species Ni_72("Ni-72", "Ni", 16, 44, 28, 72, 8520.2118, "B-", 5556.9381, 71.941785924, 2.401);
    static const Species Cu_72("Cu-72", "Cu", 14, 43, 29, 72, 8586.5256, "B-", 8362.4883, 71.935820306, 1.5);
    static const Species Zn_72("Zn-72", "Zn", 12, 42, 30, 72, 8691.8053, "B-", 442.7892, 71.926842806, 2.3);
    static const Species Ga_72("Ga-72", "Ga", 10, 41, 31, 72, 8687.0893, "B-", 3997.6263, 71.926367452, 0.878);
    static const Species Ge_72("Ge-72", "Ge", 8, 40, 32, 72, 8731.7459, "B-", -4356.1019, 71.922075824, 0.081);
    static const Species As_72("As-72", "As", 6, 39, 33, 72, 8660.3786, "B-", -361.6194, 71.926752291, 4.383);
    static const Species Se_72("Se-72", "Se", 4, 38, 34, 72, 8644.4902, "B-", -8806.4384, 71.927140506, 2.1);
    static const Species Br_72("Br-72", "Br", 2, 37, 35, 72, 8511.3126, "B-", -5121.1683, 71.936594606, 1.1);
    static const Species Kr_72("Kr-72", "Kr", 0, 36, 36, 72, 8429.3193, "B-", -15611.0, 71.942092406, 8.6);
    static const Species Rb_72("Rb-72", "Rb", -2, 35, 37, 72, 8202.0, "B-", std::numeric_limits<double>::quiet_NaN(), 71.958851, 537.0);
    static const Species Mn_73("Mn-73", "Mn", 23, 48, 25, 73, 7895.0, "B-", 17289.0, 72.992807, 644.0);
    static const Species Fe_73("Fe-73", "Fe", 21, 47, 26, 73, 8121.0, "B-", 13980.0, 72.974246, 537.0);
    static const Species Co_73("Co-73", "Co", 19, 46, 27, 73, 8302.0, "B-", 12139.0, 72.959238, 322.0);
    static const Species Ni_73("Ni-73", "Ni", 17, 45, 28, 73, 8457.6529, "B-", 8879.2856, 72.946206681, 2.601);
    static const Species Cu_73("Cu-73", "Cu", 15, 44, 29, 73, 8568.5699, "B-", 6605.9659, 72.936674376, 2.084);
    static const Species Zn_73("Zn-73", "Zn", 13, 43, 30, 73, 8648.3455, "B-", 4105.9329, 72.92958258, 2.0);
    static const Species Ga_73("Ga-73", "Ga", 11, 42, 31, 73, 8693.874, "B-", 1598.1889, 72.92517468, 1.8);
    static const Species Ge_73("Ge-73", "Ge", 9, 41, 32, 73, 8705.05, "B-", -344.7759, 72.923458954, 0.061);
    static const Species As_73("As-73", "As", 7, 40, 33, 73, 8689.6099, "B-", -2725.3604, 72.923829086, 4.136);
    static const Species Se_73("Se-73", "Se", 5, 39, 34, 73, 8641.5591, "B-", -4581.6095, 72.926754881, 7.969);
    static const Species Br_73("Br-73", "Br", 3, 38, 35, 73, 8568.0803, "B-", -7094.0287, 72.931673441, 7.237);
    static const Species Kr_73("Kr-73", "Kr", 1, 37, 36, 73, 8460.1847, "B-", -10540.1468, 72.939289193, 7.061);
    static const Species Rb_73("Rb-73", "Rb", -1, 36, 37, 73, 8305.0821, "B-", -14061.0, 72.950604506, 43.794);
    static const Species Sr_73("Sr-73", "Sr", -3, 35, 38, 73, 8102.0, "B-", std::numeric_limits<double>::quiet_NaN(), 72.9657, 430.0);
    static const Species Fe_74("Fe-74", "Fe", 22, 48, 26, 74, 8076.0, "B-", 12881.0, 73.977821, 537.0);
    static const Species Co_74("Co-74", "Co", 20, 47, 27, 74, 8239.0, "B-", 15160.0, 73.963993, 429.0);
    static const Species Ni_74("Ni-74", "Ni", 18, 46, 28, 74, 8433.0, "B-", 7306.0, 73.947718, 215.0);
    static const Species Cu_74("Cu-74", "Cu", 16, 45, 29, 74, 8521.5633, "B-", 9750.5077, 73.93987486, 6.6);
    static const Species Zn_74("Zn-74", "Zn", 14, 44, 30, 74, 8642.7547, "B-", 2292.9057, 73.92940726, 2.7);
    static const Species Ga_74("Ga-74", "Ga", 12, 43, 31, 74, 8663.1676, "B-", 5372.8249, 73.926945725, 3.214);
    static const Species Ge_74("Ge-74", "Ge", 10, 42, 32, 74, 8725.2011, "B-", -2562.3871, 73.92117776, 0.013);
    static const Species As_74("As-74", "As", 8, 41, 33, 74, 8680.002, "B-", 1353.1467, 73.923928596, 1.817);
    static const Species Se_74("Se-74", "Se", 6, 40, 34, 74, 8687.7155, "B-", -6925.0492, 73.922475933, 0.015);
    static const Species Br_74("Br-74", "Br", 4, 39, 35, 74, 8583.5615, "B-", -2956.3173, 73.929910279, 6.264);
    static const Species Kr_74("Kr-74", "Kr", 2, 38, 36, 74, 8533.039, "B-", -10415.828, 73.933084016, 2.161);
    static const Species Rb_74("Rb-74", "Rb", 0, 37, 37, 74, 8381.7123, "B-", -11089.0, 73.944265867, 3.249);
    static const Species Sr_74("Sr-74", "Sr", -2, 36, 38, 74, 8221.0, "B-", std::numeric_limits<double>::quiet_NaN(), 73.95617, 107.0);
    static const Species Fe_75("Fe-75", "Fe", 23, 49, 26, 75, 7996.0, "B-", 15861.0, 74.984219, 644.0);
    static const Species Co_75("Co-75", "Co", 21, 48, 27, 75, 8197.0, "B-", 13680.0, 74.967192, 429.0);
    static const Species Ni_75("Ni-75", "Ni", 19, 47, 28, 75, 8369.0, "B-", 10230.0, 74.952506, 215.0);
    static const Species Cu_75("Cu-75", "Cu", 17, 46, 29, 75, 8495.0801, "B-", 8088.6967, 74.941523817, 0.77);
    static const Species Zn_75("Zn-75", "Zn", 15, 45, 30, 75, 8592.4981, "B-", 5901.7231, 74.932840244, 2.1);
    static const Species Ga_75("Ga-75", "Ga", 13, 44, 31, 75, 8660.7565, "B-", 3396.3337, 74.926504484, 0.72);
    static const Species Ge_75("Ge-75", "Ge", 11, 43, 32, 75, 8695.6096, "B-", 1177.2301, 74.92285837, 0.055);
    static const Species As_75("As-75", "As", 9, 42, 33, 75, 8700.8748, "B-", -864.7139, 74.921594562, 0.948);
    static const Species Se_75("Se-75", "Se", 7, 41, 34, 75, 8678.9139, "B-", -3062.4694, 74.92252287, 0.078);
    static const Species Br_75("Br-75", "Br", 5, 40, 35, 75, 8627.6497, "B-", -4783.388, 74.925810566, 4.6);
    static const Species Kr_75("Kr-75", "Kr", 3, 39, 36, 75, 8553.4399, "B-", -7104.9299, 74.930945744, 8.7);
    static const Species Rb_75("Rb-75", "Rb", 1, 38, 37, 75, 8448.2762, "B-", -10600.0, 74.9385732, 1.266);
    static const Species Sr_75("Sr-75", "Sr", -1, 37, 38, 75, 8296.5116, "B-", -14799.0, 74.949952767, 236.183);
    static const Species Y_75("Y-75", "Y", -3, 36, 39, 75, 8089.0, "B-", std::numeric_limits<double>::quiet_NaN(), 74.96584, 322.0);
    static const Species Fe_76("Fe-76", "Fe", 24, 50, 26, 76, 7943.0, "B-", 15070.0, 75.988631, 644.0);
    static const Species Co_76("Co-76", "Co", 22, 49, 27, 76, 8131.0, "B-", 16530.0, 75.972453, 537.0);
    static const Species Ni_76("Ni-76", "Ni", 20, 48, 28, 76, 8338.0, "B-", 8791.0, 75.954707, 322.0);
    static const Species Cu_76("Cu-76", "Cu", 18, 47, 29, 76, 8443.6018, "B-", 11321.3964, 75.945268974, 0.98);
    static const Species Zn_76("Zn-76", "Zn", 16, 46, 30, 76, 8582.2735, "B-", 3993.6241, 75.933114956, 1.562);
    static const Species Ga_76("Ga-76", "Ga", 14, 45, 31, 76, 8624.5272, "B-", 6916.2501, 75.928827624, 2.1);
    static const Species Ge_76("Ge-76", "Ge", 12, 44, 32, 76, 8705.2364, "B-", -921.5145, 75.921402725, 0.019);
    static const Species As_76("As-76", "As", 10, 43, 33, 76, 8682.8172, "B-", 2960.5756, 75.922392011, 0.951);
    static const Species Se_76("Se-76", "Se", 8, 42, 34, 76, 8711.4781, "B-", -4962.881, 75.919213702, 0.017);
    static const Species Br_76("Br-76", "Br", 6, 41, 35, 76, 8635.883, "B-", -1275.3724, 75.924541574, 10.007);
    static const Species Kr_76("Kr-76", "Kr", 4, 40, 36, 76, 8608.8077, "B-", -8534.6172, 75.925910743, 4.308);
    static const Species Rb_76("Rb-76", "Rb", 2, 39, 37, 76, 8486.2161, "B-", -6231.4432, 75.935073031, 1.006);
    static const Species Sr_76("Sr-76", "Sr", 0, 38, 38, 76, 8393.9294, "B-", -15998.0, 75.94176276, 37.0);
    static const Species Y_76("Y-76", "Y", -2, 37, 39, 76, 8173.0, "B-", std::numeric_limits<double>::quiet_NaN(), 75.958937, 322.0);
    static const Species Co_77("Co-77", "Co", 23, 50, 27, 77, 8082.0, "B-", 15440.0, 76.976479, 644.0);
    static const Species Ni_77("Ni-77", "Ni", 21, 49, 28, 77, 8272.0, "B-", 11513.0, 76.959903, 429.0);
    static const Species Cu_77("Cu-77", "Cu", 19, 48, 29, 77, 8411.2501, "B-", 9926.375, 76.947543599, 1.3);
    static const Species Zn_77("Zn-77", "Zn", 17, 47, 30, 77, 8530.0037, "B-", 7203.1495, 76.936887197, 2.117);
    static const Species Ga_77("Ga-77", "Ga", 15, 46, 31, 77, 8613.3907, "B-", 5220.5176, 76.929154299, 2.6);
    static const Species Ge_77("Ge-77", "Ge", 13, 45, 32, 77, 8671.0293, "B-", 2703.4642, 76.923549843, 0.056);
    static const Species As_77("As-77", "As", 11, 44, 33, 77, 8695.9789, "B-", 683.1627, 76.920647555, 1.816);
    static const Species Se_77("Se-77", "Se", 9, 43, 34, 77, 8694.6908, "B-", -1364.6792, 76.91991415, 0.067);
    static const Species Br_77("Br-77", "Br", 7, 42, 35, 77, 8666.8073, "B-", -3065.3663, 76.921379193, 3.017);
    static const Species Kr_77("Kr-77", "Kr", 5, 41, 36, 77, 8616.837, "B-", -5338.9516, 76.924669999, 2.1);
    static const Species Rb_77("Rb-77", "Rb", 3, 40, 37, 77, 8537.3396, "B-", -7027.0566, 76.930401599, 1.4);
    static const Species Sr_77("Sr-77", "Sr", 1, 39, 38, 77, 8435.9188, "B-", -11365.0, 76.937945454, 8.5);
    static const Species Y_77("Y-77", "Y", -1, 38, 39, 77, 8278.0, "B-", -14839.0, 76.950146, 218.0);
    static const Species Zr_77("Zr-77", "Zr", -3, 37, 40, 77, 8075.0, "B-", std::numeric_limits<double>::quiet_NaN(), 76.966076, 429.0);
    static const Species Co_78("Co-78", "Co", 24, 51, 27, 78, 7997.0, "B-", 19560.0, 77.983553, 751.0);
    static const Species Ni_78("Ni-78", "Ni", 22, 50, 28, 78, 8238.0, "B-", 9910.0, 77.962555, 429.0);
    static const Species Cu_78("Cu-78", "Cu", 20, 49, 29, 78, 8354.6695, "B-", 12693.768, 77.951916524, 14.312);
    static const Species Zn_78("Zn-78", "Zn", 18, 48, 30, 78, 8507.38, "B-", 6220.8433, 77.938289204, 2.086);
    static const Species Ga_78("Ga-78", "Ga", 16, 47, 31, 78, 8577.1043, "B-", 8157.9729, 77.931610854, 1.127);
    static const Species Ge_78("Ge-78", "Ge", 14, 46, 32, 78, 8671.6636, "B-", 954.9114, 77.922852911, 3.795);
    static const Species As_78("As-78", "As", 12, 45, 33, 78, 8673.876, "B-", 4208.9819, 77.921827771, 10.498);
    static const Species Se_78("Se-78", "Se", 10, 44, 34, 78, 8717.8072, "B-", -3573.7836, 77.917309244, 0.191);
    static const Species Br_78("Br-78", "Br", 8, 43, 35, 78, 8661.9594, "B-", 726.1153, 77.921145858, 3.842);
    static const Species Kr_78("Kr-78", "Kr", 6, 42, 36, 78, 8661.2385, "B-", -7242.856, 77.920366341, 0.329);
    static const Species Rb_78("Rb-78", "Rb", 4, 41, 37, 78, 8558.3512, "B-", -3761.4779, 77.928141866, 3.475);
    static const Species Sr_78("Sr-78", "Sr", 2, 40, 38, 78, 8500.0971, "B-", -11001.0, 77.932179979, 8.0);
    static const Species Y_78("Y-78", "Y", 0, 39, 39, 78, 8349.0, "B-", -11323.0, 77.94399, 320.0);
    static const Species Zr_78("Zr-78", "Zr", -2, 38, 40, 78, 8194.0, "B-", std::numeric_limits<double>::quiet_NaN(), 77.956146, 429.0);
    static const Species Ni_79("Ni-79", "Ni", 23, 51, 28, 79, 8150.0, "B-", 14248.0, 78.969769, 537.0);
    static const Species Cu_79("Cu-79", "Cu", 21, 50, 29, 79, 8320.938, "B-", 11024.2629, 78.9544731, 112.7);
    static const Species Zn_79("Zn-79", "Zn", 19, 49, 30, 79, 8450.5825, "B-", 9116.0536, 78.942638067, 2.388);
    static const Species Ga_79("Ga-79", "Ga", 17, 48, 31, 79, 8556.0725, "B-", 6978.8242, 78.932851582, 1.296);
    static const Species Ge_79("Ge-79", "Ge", 15, 47, 32, 79, 8634.5089, "B-", 4108.9014, 78.925359506, 39.893);
    static const Species As_79("As-79", "As", 13, 46, 33, 79, 8676.6172, "B-", 2281.3849, 78.920948419, 5.716);
    static const Species Se_79("Se-79", "Se", 11, 45, 34, 79, 8695.5923, "B-", 150.6016, 78.918499252, 0.238);
    static const Species Br_79("Br-79", "Br", 9, 44, 35, 79, 8687.5956, "B-", -1625.7778, 78.918337574, 1.074);
    static const Species Kr_79("Kr-79", "Kr", 7, 43, 36, 79, 8657.113, "B-", -3639.5114, 78.920082919, 3.736);
    static const Species Rb_79("Rb-79", "Rb", 5, 42, 37, 79, 8601.1401, "B-", -5323.114, 78.923990095, 2.085);
    static const Species Sr_79("Sr-79", "Sr", 3, 41, 38, 79, 8523.8558, "B-", -7676.7291, 78.929704692, 7.967);
    static const Species Y_79("Y-79", "Y", 1, 40, 39, 79, 8416.7788, "B-", -11033.0, 78.937946, 86.0);
    static const Species Zr_79("Zr-79", "Zr", -1, 39, 40, 79, 8267.0, "B-", -15120.0, 78.94979, 322.0);
    static const Species Nb_79("Nb-79", "Nb", -3, 38, 41, 79, 8066.0, "B-", std::numeric_limits<double>::quiet_NaN(), 78.966022, 537.0);
    static const Species Ni_80("Ni-80", "Ni", 24, 52, 28, 80, 8088.0, "B-", 13440.0, 79.975051, 644.0);
    static const Species Cu_80("Cu-80", "Cu", 22, 51, 29, 80, 8246.0, "B-", 14969.0, 79.960623, 322.0);
    static const Species Zn_80("Zn-80", "Zn", 20, 50, 30, 80, 8423.5457, "B-", 7575.0553, 79.944552929, 2.774);
    static const Species Ga_80("Ga-80", "Ga", 18, 49, 31, 80, 8508.4545, "B-", 10311.6397, 79.936420773, 3.103);
    static const Species Ge_80("Ge-80", "Ge", 16, 48, 32, 80, 8627.5707, "B-", 2679.2869, 79.925350773, 2.205);
    static const Species As_80("As-80", "As", 14, 47, 33, 80, 8651.2824, "B-", 5544.8861, 79.92247444, 3.578);
    static const Species Se_80("Se-80", "Se", 12, 46, 34, 80, 8710.8142, "B-", -1870.4623, 79.916521761, 1.016);
    static const Species Br_80("Br-80", "Br", 10, 45, 35, 80, 8677.6541, "B-", 2004.4299, 79.918529784, 1.065);
    static const Species Kr_80("Kr-80", "Kr", 8, 44, 36, 80, 8692.9301, "B-", -5717.9785, 79.91637794, 0.745);
    static const Species Rb_80("Rb-80", "Rb", 6, 43, 37, 80, 8611.676, "B-", -1864.009, 79.922516442, 2.0);
    static const Species Sr_80("Sr-80", "Sr", 4, 42, 38, 80, 8578.5966, "B-", -9163.305, 79.924517538, 3.718);
    static const Species Y_80("Y-80", "Y", 2, 41, 39, 80, 8454.2759, "B-", -6388.0, 79.93435475, 6.701);
    static const Species Zr_80("Zr-80", "Zr", 0, 40, 40, 80, 8365.0, "B-", -16339.0, 79.941213, 322.0);
    static const Species Nb_80("Nb-80", "Nb", -2, 39, 41, 80, 8151.0, "B-", std::numeric_limits<double>::quiet_NaN(), 79.958754, 429.0);
    static const Species Ni_81("Ni-81", "Ni", 25, 53, 28, 81, 8000.0, "B-", 15820.0, 80.982727, 751.0);
    static const Species Cu_81("Cu-81", "Cu", 23, 52, 29, 81, 8185.0, "B-", 14289.0, 80.965743, 322.0);
    static const Species Zn_81("Zn-81", "Zn", 21, 51, 30, 81, 8351.9262, "B-", 11428.2924, 80.950402617, 5.4);
    static const Species Ga_81("Ga-81", "Ga", 19, 50, 31, 81, 8483.3576, "B-", 8663.7335, 80.938133841, 3.503);
    static const Species Ge_81("Ge-81", "Ge", 17, 49, 32, 81, 8580.6587, "B-", 6241.6189, 80.928832941, 2.205);
    static const Species As_81("As-81", "As", 15, 48, 33, 81, 8648.0571, "B-", 3855.705, 80.922132288, 2.838);
    static const Species Se_81("Se-81", "Se", 13, 47, 34, 81, 8685.9998, "B-", 1588.0317, 80.917993019, 1.049);
    static const Species Br_81("Br-81", "Br", 11, 46, 35, 81, 8695.9465, "B-", -280.8517, 80.916288197, 1.049);
    static const Species Kr_81("Kr-81", "Kr", 9, 45, 36, 81, 8682.8206, "B-", -2239.4954, 80.916589703, 1.152);
    static const Species Rb_81("Rb-81", "Rb", 7, 44, 37, 81, 8645.5139, "B-", -3928.5695, 80.9189939, 5.265);
    static const Species Sr_81("Sr-81", "Sr", 5, 43, 38, 81, 8587.3545, "B-", -5815.2156, 80.923211393, 3.358);
    static const Species Y_81("Y-81", "Y", 3, 42, 39, 81, 8505.9031, "B-", -8188.5003, 80.929454283, 5.802);
    static const Species Zr_81("Zr-81", "Zr", 1, 41, 40, 81, 8395.1519, "B-", -11164.0, 80.938245, 99.0);
    static const Species Nb_81("Nb-81", "Nb", -1, 40, 41, 81, 8248.0, "B-", -14900.0, 80.95023, 429.0);
    static const Species Mo_81("Mo-81", "Mo", -3, 39, 42, 81, 8054.0, "B-", std::numeric_limits<double>::quiet_NaN(), 80.966226, 537.0);
    static const Species Ni_82("Ni-82", "Ni", 26, 54, 28, 82, 7935.0, "B-", 15010.0, 81.988492, 859.0);
    static const Species Cu_82("Cu-82", "Cu", 24, 53, 29, 82, 8108.0, "B-", 16584.0, 81.972378, 429.0);
    static const Species Zn_82("Zn-82", "Zn", 22, 52, 30, 82, 8301.1175, "B-", 10616.7652, 81.954574097, 3.3);
    static const Species Ga_82("Ga-82", "Ga", 20, 51, 31, 82, 8421.0494, "B-", 12484.3497, 81.943176531, 2.604);
    static const Species Ge_82("Ge-82", "Ge", 18, 50, 32, 82, 8563.7567, "B-", 4690.3523, 81.929774031, 2.405);
    static const Species As_82("As-82", "As", 16, 49, 33, 82, 8611.4153, "B-", 7488.4677, 81.924738731, 4.003);
    static const Species Se_82("Se-82", "Se", 14, 48, 34, 82, 8693.1973, "B-", -95.2184, 81.916699531, 0.5);
    static const Species Br_82("Br-82", "Br", 12, 47, 35, 82, 8682.4953, "B-", 3093.1185, 81.916801752, 1.042);
    static const Species Kr_82("Kr-82", "Kr", 10, 46, 36, 82, 8710.6754, "B-", -4403.9824, 81.91348115368, 0.00591);
    static const Species Rb_82("Rb-82", "Rb", 8, 45, 37, 82, 8647.4275, "B-", -177.7503, 81.918209023, 3.23);
    static const Species Sr_82("Sr-82", "Sr", 6, 44, 38, 82, 8635.719, "B-", -7945.965, 81.918399845, 6.432);
    static const Species Y_82("Y-82", "Y", 4, 43, 39, 82, 8529.2762, "B-", -4450.0341, 81.926930189, 5.902);
    static const Species Zr_82("Zr-82", "Zr", 2, 42, 40, 82, 8465.4666, "B-", -11804.0, 81.931707497, 1.7);
    static const Species Nb_82("Nb-82", "Nb", 0, 41, 41, 82, 8312.0, "B-", -11440.0, 81.94438, 322.0);
    static const Species Mo_82("Mo-82", "Mo", -2, 40, 42, 82, 8163.0, "B-", std::numeric_limits<double>::quiet_NaN(), 81.956661, 429.0);
    static const Species Cu_83("Cu-83", "Cu", 25, 54, 29, 83, 8044.0, "B-", 15900.0, 82.97811, 537.0);
    static const Species Zn_83("Zn-83", "Zn", 23, 53, 30, 83, 8226.0, "B-", 12967.0, 82.961041, 322.0);
    static const Species Ga_83("Ga-83", "Ga", 21, 52, 31, 83, 8372.5756, "B-", 11719.3136, 82.9471203, 2.804);
    static const Species Ge_83("Ge-83", "Ge", 19, 51, 32, 83, 8504.3462, "B-", 8692.8893, 82.9345391, 2.604);
    static const Species As_83("As-83", "As", 17, 50, 33, 83, 8599.654, "B-", 5671.2117, 82.9252069, 3.004);
    static const Species Se_83("Se-83", "Se", 15, 49, 34, 83, 8658.556, "B-", 3673.178, 82.919118604, 3.259);
    static const Species Br_83("Br-83", "Br", 13, 48, 35, 83, 8693.3852, "B-", 976.9222, 82.915175285, 4.073);
    static const Species Kr_83("Kr-83", "Kr", 11, 47, 36, 83, 8695.7295, "B-", -920.0039, 82.914126516, 0.009);
    static const Species Rb_83("Rb-83", "Rb", 9, 46, 37, 83, 8675.2193, "B-", -2273.0239, 82.915114181, 2.5);
    static const Species Sr_83("Sr-83", "Sr", 7, 45, 38, 83, 8638.4076, "B-", -4591.9435, 82.917554372, 7.336);
    static const Species Y_83("Y-83", "Y", 5, 44, 39, 83, 8573.6571, "B-", -6294.0125, 82.922484026, 20.0);
    static const Species Zr_83("Zr-83", "Zr", 3, 43, 40, 83, 8488.3997, "B-", -8298.7493, 82.929240926, 6.902);
    static const Species Nb_83("Nb-83", "Nb", 1, 42, 41, 83, 8378.9889, "B-", -11273.0, 82.93815, 174.0);
    static const Species Mo_83("Mo-83", "Mo", -1, 41, 42, 83, 8234.0, "B-", -15020.0, 82.950252, 430.0);
    static const Species Tc_83("Tc-83", "Tc", -3, 40, 43, 83, 8043.0, "B-", std::numeric_limits<double>::quiet_NaN(), 82.966377, 537.0);
    static const Species Cu_84("Cu-84", "Cu", 26, 55, 29, 84, 7965.0, "B-", 18110.0, 83.985271, 537.0);
    static const Species Zn_84("Zn-84", "Zn", 24, 54, 30, 84, 8171.0, "B-", 12264.0, 83.965829, 429.0);
    static const Species Ga_84("Ga-84", "Ga", 22, 53, 31, 84, 8307.525, "B-", 14054.2989, 83.952663, 32.0);
    static const Species Ge_84("Ge-84", "Ge", 20, 52, 32, 84, 8465.5244, "B-", 7705.1329, 83.93757509, 3.403);
    static const Species As_84("As-84", "As", 18, 51, 33, 84, 8547.9385, "B-", 10094.1624, 83.92930329000001, 3.403);
    static const Species Se_84("Se-84", "Se", 16, 50, 34, 84, 8658.7935, "B-", 1835.3638, 83.918466761, 2.105);
    static const Species Br_84("Br-84", "Br", 14, 49, 35, 84, 8671.3294, "B-", 4656.251, 83.916496417, 27.622);
    static const Species Kr_84("Kr-84", "Kr", 12, 48, 36, 84, 8717.4473, "B-", -2680.3708, 83.91149772708, 0.0041);
    static const Species Rb_84("Rb-84", "Rb", 10, 47, 37, 84, 8676.2244, "B-", 890.6058, 83.914375223, 2.355);
    static const Species Sr_84("Sr-84", "Sr", 8, 46, 38, 84, 8677.5132, "B-", -6755.1411, 83.913419118, 1.334);
    static const Species Y_84("Y-84", "Y", 6, 45, 39, 84, 8587.7812, "B-", -2472.7471, 83.92067106, 4.615);
    static const Species Zr_84("Zr-84", "Zr", 4, 44, 40, 84, 8549.0301, "B-", -10227.8497, 83.923325663, 5.903);
    static const Species Nb_84("Nb-84", "Nb", 2, 43, 41, 84, 8417.9563, "B-", -7024.0, 83.934305711, 0.43);
    static const Species Mo_84("Mo-84", "Mo", 0, 42, 42, 84, 8325.0, "B-", -16470.0, 83.941846, 320.0);
    static const Species Tc_84("Tc-84", "Tc", -2, 41, 43, 84, 8120.0, "B-", std::numeric_limits<double>::quiet_NaN(), 83.959527, 429.0);
    static const Species Zn_85("Zn-85", "Zn", 25, 55, 30, 85, 8090.0, "B-", 14644.0, 84.973054, 537.0);
    static const Species Ga_85("Ga-85", "Ga", 23, 54, 31, 85, 8253.5687, "B-", 13379.3679, 84.957333, 40.0);
    static const Species Ge_85("Ge-85", "Ge", 21, 53, 32, 85, 8401.7689, "B-", 10065.7253, 84.942969658, 4.003);
    static const Species As_85("As-85", "As", 19, 52, 33, 85, 8510.9851, "B-", 9224.4929, 84.932163658, 3.304);
    static const Species Se_85("Se-85", "Se", 17, 51, 34, 85, 8610.3045, "B-", 6161.8335, 84.922260758, 2.804);
    static const Species Br_85("Br-85", "Br", 15, 50, 35, 85, 8673.5926, "B-", 2904.8622, 84.915645758, 3.304);
    static const Species Kr_85("Kr-85", "Kr", 13, 49, 36, 85, 8698.5633, "B-", 687.0, 84.91252726, 2.147);
    static const Species Rb_85("Rb-85", "Rb", 11, 48, 37, 85, 8697.4416, "B-", -1064.051, 84.91178973604, 0.00537);
    static const Species Sr_85("Sr-85", "Sr", 9, 47, 38, 85, 8675.7193, "B-", -3261.1584, 84.912932041, 3.02);
    static const Species Y_85("Y-85", "Y", 7, 46, 39, 85, 8628.1486, "B-", -4666.9352, 84.916433039, 20.36);
    static const Species Zr_85("Zr-85", "Zr", 5, 45, 40, 85, 8564.0394, "B-", -6895.512, 84.921443199, 6.902);
    static const Species Nb_85("Nb-85", "Nb", 3, 44, 41, 85, 8473.7117, "B-", -8769.9238, 84.928845836, 4.4);
    static const Species Mo_85("Mo-85", "Mo", 1, 43, 42, 85, 8361.332, "B-", -11660.0, 84.938260736, 17.0);
    static const Species Tc_85("Tc-85", "Tc", -1, 42, 43, 85, 8215.0, "B-", -15220.0, 84.950778, 429.0);
    static const Species Ru_85("Ru-85", "Ru", -3, 41, 44, 85, 8027.0, "B-", std::numeric_limits<double>::quiet_NaN(), 84.967117, 537.0);
    static const Species Zn_86("Zn-86", "Zn", 26, 56, 30, 86, 8032.0, "B-", 13699.0, 85.978463, 537.0);
    static const Species Ga_86("Ga-86", "Ga", 24, 55, 31, 86, 8182.0, "B-", 15640.0, 85.963757, 429.0);
    static const Species Ge_86("Ge-86", "Ge", 22, 54, 32, 86, 8354.63, "B-", 9562.2229, 85.946967, 470.0);
    static const Species As_86("As-86", "As", 20, 53, 33, 86, 8456.7215, "B-", 11541.0256, 85.936701532, 3.703);
    static const Species Se_86("Se-86", "Se", 18, 52, 34, 86, 8581.8224, "B-", 5129.086, 85.924311732, 2.705);
    static const Species Br_86("Br-86", "Br", 16, 51, 35, 86, 8632.3659, "B-", 7633.4147, 85.918805432, 3.304);
    static const Species Kr_86("Kr-86", "Kr", 14, 50, 36, 86, 8712.0295, "B-", -518.6734, 85.91061062468, 0.00399);
    static const Species Rb_86("Rb-86", "Rb", 12, 49, 37, 86, 8696.9014, "B-", 1776.0972, 85.911167443, 0.214);
    static const Species Sr_86("Sr-86", "Sr", 10, 48, 38, 86, 8708.4566, "B-", -5240.0, 85.90926072473, 0.00563);
    static const Species Y_86("Y-86", "Y", 8, 47, 39, 86, 8638.4293, "B-", -1314.0763, 85.914886095, 15.182);
    static const Species Zr_86("Zr-86", "Zr", 6, 46, 40, 86, 8614.0523, "B-", -8834.9627, 85.916296814, 3.827);
    static const Species Nb_86("Nb-86", "Nb", 4, 45, 41, 86, 8502.2231, "B-", -5023.1337, 85.925781536, 5.903);
    static const Species Mo_86("Mo-86", "Mo", 2, 44, 42, 86, 8434.7175, "B-", -12541.0, 85.931174092, 3.147);
    static const Species Tc_86("Tc-86", "Tc", 0, 43, 43, 86, 8280.0, "B-", -11800.0, 85.944637, 322.0);
    static const Species Ru_86("Ru-86", "Ru", -2, 42, 44, 86, 8133.0, "B-", std::numeric_limits<double>::quiet_NaN(), 85.957305, 429.0);
    static const Species Ga_87("Ga-87", "Ga", 25, 56, 31, 87, 8124.0, "B-", 14720.0, 86.969007, 537.0);
    static const Species Ge_87("Ge-87", "Ge", 23, 55, 32, 87, 8285.0, "B-", 12028.0, 86.953204, 322.0);
    static const Species As_87("As-87", "As", 21, 54, 33, 87, 8413.8521, "B-", 10808.2192, 86.940291716, 3.204);
    static const Species Se_87("Se-87", "Se", 19, 53, 34, 87, 8529.092, "B-", 7465.5526, 86.928688616, 2.405);
    static const Species Br_87("Br-87", "Br", 17, 52, 35, 87, 8605.9105, "B-", 6817.8455, 86.920674016, 3.404);
    static const Species Kr_87("Kr-87", "Kr", 15, 51, 36, 87, 8675.284, "B-", 3888.2706, 86.913354759, 0.264);
    static const Species Rb_87("Rb-87", "Rb", 13, 50, 37, 87, 8710.9843, "B-", 282.2749, 86.909180529, 0.006);
    static const Species Sr_87("Sr-87", "Sr", 11, 49, 38, 87, 8705.2363, "B-", -1861.6894, 86.90887749454, 0.0055);
    static const Species Y_87("Y-87", "Y", 9, 48, 39, 87, 8674.8451, "B-", -3671.2405, 86.9108761, 1.21);
    static const Species Zr_87("Zr-87", "Zr", 7, 47, 40, 87, 8623.6545, "B-", -5472.6536, 86.914817338, 4.45);
    static const Species Nb_87("Nb-87", "Nb", 5, 46, 41, 87, 8551.7579, "B-", -6989.6757, 86.920692473, 7.302);
    static const Species Mo_87("Mo-87", "Mo", 3, 45, 42, 87, 8462.4243, "B-", -9194.7656, 86.928196198, 3.067);
    static const Species Tc_87("Tc-87", "Tc", 1, 44, 43, 87, 8347.7449, "B-", -11960.0, 86.938067185, 4.5);
    static const Species Ru_87("Ru-87", "Ru", -1, 43, 44, 87, 8201.0, "B-", std::numeric_limits<double>::quiet_NaN(), 86.950907, 429.0);
    static const Species Ga_88("Ga-88", "Ga", 26, 57, 31, 88, 8050.0, "B-", 17129.0, 87.975963, 537.0);
    static const Species Ge_88("Ge-88", "Ge", 24, 56, 32, 88, 8236.0, "B-", 10930.0, 87.957574, 429.0);
    static const Species As_88("As-88", "As", 22, 55, 33, 88, 8351.0, "B-", 13434.0, 87.94584, 215.0);
    static const Species Se_88("Se-88", "Se", 20, 54, 34, 88, 8495.0045, "B-", 6831.764, 87.93141749, 3.604);
    static const Species Br_88("Br-88", "Br", 18, 53, 35, 88, 8563.7479, "B-", 8975.3282, 87.92408329, 3.404);
    static const Species Kr_88("Kr-88", "Kr", 16, 52, 36, 88, 8656.8499, "B-", 2917.709, 87.914447879, 2.8);
    static const Species Rb_88("Rb-88", "Rb", 14, 51, 37, 88, 8681.1154, "B-", 5312.6243, 87.91131559, 0.17);
    static const Species Sr_88("Sr-88", "Sr", 12, 50, 38, 88, 8732.5958, "B-", -3622.6, 87.905612253, 0.006);
    static const Species Y_88("Y-88", "Y", 10, 49, 39, 88, 8682.5396, "B-", -670.1549, 87.909501274, 1.61);
    static const Species Zr_88("Zr-88", "Zr", 8, 48, 40, 88, 8666.0339, "B-", -7457.3187, 87.910220715, 5.8);
    static const Species Nb_88("Nb-88", "Nb", 6, 47, 41, 88, 8572.4013, "B-", -3485.0021, 87.918226476, 62.059);
    static const Species Mo_88("Mo-88", "Mo", 4, 46, 42, 88, 8523.9087, "B-", -11016.2515, 87.921967779, 4.1);
    static const Species Tc_88("Tc-88", "Tc", 2, 45, 43, 88, 8389.8338, "B-", -7331.0, 87.933794211, 4.4);
    static const Species Ru_88("Ru-88", "Ru", 0, 44, 44, 88, 8298.0, "B-", -17479.0, 87.941664, 322.0);
    static const Species Rh_88("Rh-88", "Rh", -2, 43, 45, 88, 8090.0, "B-", std::numeric_limits<double>::quiet_NaN(), 87.960429, 429.0);
    static const Species Ge_89("Ge-89", "Ge", 25, 57, 32, 89, 8161.0, "B-", 13490.0, 88.96453, 429.0);
    static const Species As_89("As-89", "As", 23, 56, 33, 89, 8304.0, "B-", 12462.0, 88.950048, 322.0);
    static const Species Se_89("Se-89", "Se", 21, 55, 34, 89, 8435.2799, "B-", 9281.873, 88.936669058, 4.003);
    static const Species Br_89("Br-89", "Br", 19, 54, 35, 89, 8530.7802, "B-", 8261.5231, 88.926704558, 3.504);
    static const Species Kr_89("Kr-89", "Kr", 17, 53, 36, 89, 8614.8158, "B-", 5176.6042, 88.917835449, 2.3);
    static const Species Rb_89("Rb-89", "Rb", 15, 52, 37, 89, 8664.1895, "B-", 4496.6278, 88.912278136, 5.825);
    static const Species Sr_89("Sr-89", "Sr", 13, 51, 38, 89, 8705.923, "B-", 1502.1757, 88.907450808, 0.098);
    static const Species Y_89("Y-89", "Y", 11, 50, 39, 89, 8714.011, "B-", -2833.2285, 88.905838156, 0.363);
    static const Species Zr_89("Zr-89", "Zr", 9, 49, 40, 89, 8673.3865, "B-", -4252.2191, 88.908879751, 2.983);
    static const Species Nb_89("Nb-89", "Nb", 7, 48, 41, 89, 8616.8184, "B-", -5610.8105, 88.913444696, 25.367);
    static const Species Mo_89("Mo-89", "Mo", 5, 47, 42, 89, 8544.9851, "B-", -7620.0875, 88.919468149, 4.2);
    static const Species Tc_89("Tc-89", "Tc", 3, 46, 43, 89, 8450.5758, "B-", -9025.4327, 88.927648649, 4.1);
    static const Species Ru_89("Ru-89", "Ru", 1, 45, 44, 89, 8340.376, "B-", -12719.0, 88.937337849, 26.0);
    static const Species Rh_89("Rh-89", "Rh", -1, 44, 45, 89, 8189.0, "B-", std::numeric_limits<double>::quiet_NaN(), 88.950992, 387.0);
    static const Species Ge_90("Ge-90", "Ge", 26, 58, 32, 90, 8109.0, "B-", 12520.0, 89.969436, 537.0);
    static const Species As_90("As-90", "As", 24, 57, 33, 90, 8240.0, "B-", 14810.0, 89.955995, 429.0);
    static const Species Se_90("Se-90", "Se", 22, 56, 34, 90, 8395.7672, "B-", 8200.0834, 89.940096, 354.0);
    static const Species Br_90("Br-90", "Br", 20, 55, 35, 90, 8478.1865, "B-", 10958.9533, 89.931292848, 3.604);
    static const Species Kr_90("Kr-90", "Kr", 18, 54, 36, 90, 8591.2599, "B-", 4406.3133, 89.919527929, 2.0);
    static const Species Rb_90("Rb-90", "Rb", 16, 53, 37, 90, 8631.5262, "B-", 6585.3721, 89.914797557, 6.926);
    static const Species Sr_90("Sr-90", "Sr", 14, 52, 38, 90, 8696.0043, "B-", 545.9674, 89.90772787, 1.555);
    static const Species Y_90("Y-90", "Y", 12, 51, 39, 90, 8693.3778, "B-", 2275.635, 89.907141749, 0.379);
    static const Species Zr_90("Zr-90", "Zr", 10, 50, 40, 90, 8709.9699, "B-", -6111.0165, 89.904698755, 0.126);
    static const Species Nb_90("Nb-90", "Nb", 8, 49, 41, 90, 8633.377, "B-", -2489.0165, 89.91125920100001, 3.561);
    static const Species Mo_90("Mo-90", "Mo", 6, 48, 42, 90, 8597.0285, "B-", -9447.8181, 89.91393127, 3.717);
    static const Species Tc_90("Tc-90", "Tc", 4, 47, 43, 90, 8483.36, "B-", -5840.8951, 89.924073919, 1.1);
    static const Species Ru_90("Ru-90", "Ru", 2, 46, 44, 90, 8409.7684, "B-", -13250.0, 89.930344378, 4.004);
    static const Species Rh_90("Rh-90", "Rh", 0, 45, 45, 90, 8254.0, "B-", -11924.0, 89.944569, 215.0);
    static const Species Pd_90("Pd-90", "Pd", -2, 44, 46, 90, 8113.0, "B-", std::numeric_limits<double>::quiet_NaN(), 89.95737, 429.0);
    static const Species As_91("As-91", "As", 25, 58, 33, 91, 8189.0, "B-", 14080.0, 90.960816, 429.0);
    static const Species Se_91("Se-91", "Se", 23, 57, 34, 91, 8334.8382, "B-", 10527.1716, 90.9457, 465.0);
    static const Species Br_91("Br-91", "Br", 21, 56, 35, 91, 8441.9242, "B-", 9866.6724, 90.934398617, 3.804);
    static const Species Kr_91("Kr-91", "Kr", 19, 55, 36, 91, 8541.7519, "B-", 6771.0748, 90.923806309, 2.4);
    static const Species Rb_91("Rb-91", "Rb", 17, 54, 37, 91, 8607.5621, "B-", 5906.901, 90.916537261, 8.375);
    static const Species Sr_91("Sr-91", "Sr", 15, 53, 38, 91, 8663.8759, "B-", 2699.3714, 90.910195942, 5.853);
    static const Species Y_91("Y-91", "Y", 13, 52, 39, 91, 8684.9421, "B-", 1544.271, 90.907298048, 1.978);
    static const Species Zr_91("Zr-91", "Zr", 11, 51, 40, 91, 8693.3149, "B-", -1257.5644, 90.905640205, 0.101);
    static const Species Nb_91("Nb-91", "Nb", 9, 50, 41, 91, 8670.8983, "B-", -4429.1934, 90.906990256, 3.14);
    static const Species Mo_91("Mo-91", "Mo", 7, 49, 42, 91, 8613.6286, "B-", -6222.1768, 90.91174519, 6.696);
    static const Species Tc_91("Tc-91", "Tc", 5, 48, 43, 91, 8536.6558, "B-", -7746.8246, 90.918424972, 2.536);
    static const Species Ru_91("Ru-91", "Ru", 3, 47, 44, 91, 8442.9287, "B-", -9670.0, 90.92674153, 2.384);
    static const Species Rh_91("Rh-91", "Rh", 1, 46, 45, 91, 8328.0, "B-", -12400.0, 90.937123, 320.0);
    static const Species Pd_91("Pd-91", "Pd", -1, 45, 46, 91, 8183.0, "B-", std::numeric_limits<double>::quiet_NaN(), 90.950435, 454.0);
    static const Species As_92("As-92", "As", 26, 59, 33, 92, 8121.0, "B-", 16344.0, 91.967386, 537.0);
    static const Species Se_92("Se-92", "Se", 24, 58, 34, 92, 8290.0, "B-", 9509.0, 91.94984, 429.0);
    static const Species Br_92("Br-92", "Br", 22, 57, 35, 92, 8384.9123, "B-", 12536.5161, 91.939631595, 7.202);
    static const Species Kr_92("Kr-92", "Kr", 20, 56, 36, 92, 8512.675, "B-", 6003.121, 91.926173092, 2.9);
    static const Species Rb_92("Rb-92", "Rb", 18, 55, 37, 92, 8569.4225, "B-", 8094.9212, 91.919728477, 6.573);
    static const Species Sr_92("Sr-92", "Sr", 16, 54, 38, 92, 8648.907, "B-", 1949.1237, 91.911038222, 3.675);
    static const Species Y_92("Y-92", "Y", 14, 53, 39, 92, 8661.5894, "B-", 3642.5294, 91.908945752, 9.798);
    static const Species Zr_92("Zr-92", "Zr", 12, 52, 40, 92, 8692.6783, "B-", -2005.7335, 91.905035336, 0.101);
    static const Species Nb_92("Nb-92", "Nb", 10, 51, 41, 92, 8662.3731, "B-", 355.2968, 91.90718858, 1.915);
    static const Species Mo_92("Mo-92", "Mo", 8, 50, 42, 92, 8657.7312, "B-", -7882.8841, 91.906807153, 0.168);
    static const Species Tc_92("Tc-92", "Tc", 6, 49, 43, 92, 8563.544, "B-", -4624.4922, 91.915269777, 3.33);
    static const Species Ru_92("Ru-92", "Ru", 4, 48, 44, 92, 8504.774, "B-", -11302.1155, 91.920234373, 2.917);
    static const Species Rh_92("Rh-92", "Rh", 2, 47, 45, 92, 8373.4211, "B-", -8220.0, 91.932367692, 4.7);
    static const Species Pd_92("Pd-92", "Pd", 0, 46, 46, 92, 8275.5695, "B-", -17249.0, 91.941192225, 370.402);
    static const Species Ag_92("Ag-92", "Ag", -2, 45, 47, 92, 8080.0, "B-", std::numeric_limits<double>::quiet_NaN(), 91.95971, 429.0);
    static const Species Se_93("Se-93", "Se", 25, 59, 34, 93, 8225.0, "B-", 12030.0, 92.956135, 429.0);
    static const Species Br_93("Br-93", "Br", 23, 58, 35, 93, 8345.5986, "B-", 11245.7673, 92.94322, 462.5);
    static const Species Kr_93("Kr-93", "Kr", 21, 57, 36, 93, 8458.1085, "B-", 8483.8977, 92.931147172, 2.7);
    static const Species Rb_93("Rb-93", "Rb", 19, 56, 37, 93, 8540.9209, "B-", 7465.9434, 92.922039334, 8.406);
    static const Species Sr_93("Sr-93", "Sr", 17, 55, 38, 93, 8612.7875, "B-", 4141.3118, 92.914024314, 8.109);
    static const Species Y_93("Y-93", "Y", 15, 54, 39, 93, 8648.9054, "B-", 2894.8723, 92.909578434, 11.259);
    static const Species Zr_93("Zr-93", "Zr", 13, 53, 40, 93, 8671.6207, "B-", 90.8123, 92.906470661, 0.489);
    static const Species Nb_93("Nb-93", "Nb", 11, 52, 41, 93, 8664.1849, "B-", -405.7609, 92.90637317, 1.599);
    static const Species Mo_93("Mo-93", "Mo", 9, 51, 42, 93, 8651.4095, "B-", -3200.9629, 92.906808772, 0.193);
    static const Species Tc_93("Tc-93", "Tc", 7, 50, 43, 93, 8608.5782, "B-", -6389.3929, 92.910245147, 1.086);
    static const Species Ru_93("Ru-93", "Ru", 5, 49, 44, 93, 8531.4627, "B-", -8204.9136, 92.917104442, 2.216);
    static const Species Rh_93("Rh-93", "Rh", 3, 48, 45, 93, 8434.8255, "B-", -10030.0, 92.925912778, 2.821);
    static const Species Pd_93("Pd-93", "Pd", 1, 47, 46, 93, 8318.5637, "B-", -12582.0, 92.936680426, 397.221);
    static const Species Ag_93("Ag-93", "Ag", -1, 46, 47, 93, 8175.0, "B-", std::numeric_limits<double>::quiet_NaN(), 92.950188, 430.0);
    static const Species Se_94("Se-94", "Se", 26, 60, 34, 94, 8180.0, "B-", 10846.0, 93.96049, 537.0);
    static const Species Br_94("Br-94", "Br", 24, 59, 35, 94, 8287.0, "B-", 13698.0, 93.948846, 215.0);
    static const Species Kr_94("Kr-94", "Kr", 22, 58, 36, 94, 8424.3318, "B-", 7215.0114, 93.934140452, 13.0);
    static const Species Rb_94("Rb-94", "Rb", 20, 57, 37, 94, 8492.7644, "B-", 10282.9297, 93.926394819, 2.177);
    static const Species Sr_94("Sr-94", "Sr", 18, 56, 38, 94, 8593.8344, "B-", 3505.7517, 93.915355641, 1.785);
    static const Species Y_94("Y-94", "Y", 16, 55, 39, 94, 8622.8068, "B-", 4917.8589, 93.911592062, 6.849);
    static const Species Zr_94("Zr-94", "Zr", 14, 54, 40, 94, 8666.8016, "B-", -900.2684, 93.906312523, 0.175);
    static const Species Nb_94("Nb-94", "Nb", 12, 53, 41, 94, 8648.9014, "B-", 2045.0163, 93.907279001, 1.6);
    static const Species Mo_94("Mo-94", "Mo", 10, 52, 42, 94, 8662.3341, "B-", -4255.7476, 93.905083586, 0.151);
    static const Species Tc_94("Tc-94", "Tc", 8, 51, 43, 94, 8608.7373, "B-", -1574.7296, 93.909652319, 4.37);
    static const Species Ru_94("Ru-94", "Ru", 6, 50, 44, 94, 8583.662, "B-", -9675.9789, 93.91134286, 3.374);
    static const Species Rh_94("Rh-94", "Rh", 4, 49, 45, 94, 8472.4033, "B-", -6805.3428, 93.92173045, 3.627);
    static const Species Pd_94("Pd-94", "Pd", 2, 48, 46, 94, 8391.6832, "B-", -13700.0, 93.929036286, 4.602);
    static const Species Ag_94("Ag-94", "Ag", 0, 47, 47, 94, 8238.0, "B-", -11962.0, 93.943744, 429.0);
    static const Species Cd_94("Cd-94", "Cd", -2, 46, 48, 94, 8102.0, "B-", std::numeric_limits<double>::quiet_NaN(), 93.956586, 537.0);
    static const Species Se_95("Se-95", "Se", 27, 61, 34, 95, 8112.0, "B-", 13390.0, 94.9673, 537.0);
    static const Species Br_95("Br-95", "Br", 25, 60, 35, 95, 8245.0, "B-", 12309.0, 94.952925, 322.0);
    static const Species Kr_95("Kr-95", "Kr", 23, 59, 36, 95, 8365.9963, "B-", 9731.3868, 94.939710922, 20.0);
    static const Species Rb_95("Rb-95", "Rb", 21, 58, 37, 95, 8460.1967, "B-", 9226.9772, 94.929263849, 21.733);
    static const Species Sr_95("Sr-95", "Sr", 19, 57, 38, 95, 8549.0875, "B-", 6090.6528, 94.919358282, 6.237);
    static const Species Y_95("Y-95", "Y", 17, 56, 39, 95, 8604.9644, "B-", 4452.0031, 94.912819697, 7.277);
    static const Species Zr_95("Zr-95", "Zr", 15, 55, 40, 95, 8643.5924, "B-", 1126.3312, 94.908040276, 0.933);
    static const Species Nb_95("Nb-95", "Nb", 13, 54, 41, 95, 8647.2133, "B-", 925.6009, 94.90683111, 0.545);
    static const Species Mo_95("Mo-95", "Mo", 11, 53, 42, 95, 8648.7212, "B-", -1690.5175, 94.905837436, 0.132);
    static const Species Tc_95("Tc-95", "Tc", 9, 52, 43, 95, 8622.6911, "B-", -2563.5961, 94.907652281, 5.453);
    static const Species Ru_95("Ru-95", "Ru", 7, 51, 44, 95, 8587.4706, "B-", -5117.1423, 94.910404415, 10.2);
    static const Species Rh_95("Rh-95", "Rh", 5, 50, 45, 95, 8525.3707, "B-", -8374.7035, 94.915897893, 4.171);
    static const Species Pd_95("Pd-95", "Pd", 3, 49, 46, 95, 8428.9807, "B-", -10060.0, 94.924888506, 3.253);
    static const Species Ag_95("Ag-95", "Ag", 1, 48, 47, 95, 8315.0, "B-", -12850.0, 94.935688, 429.0);
    static const Species Cd_95("Cd-95", "Cd", -1, 47, 48, 95, 8171.0, "B-", std::numeric_limits<double>::quiet_NaN(), 94.949483, 607.0);
    static const Species Br_96("Br-96", "Br", 26, 61, 35, 96, 8184.0, "B-", 14872.0, 95.95898, 322.0);
    static const Species Kr_96("Kr-96", "Kr", 24, 60, 36, 96, 8330.8721, "B-", 8272.6693, 95.943014473, 20.695);
    static const Species Rb_96("Rb-96", "Rb", 22, 59, 37, 96, 8408.8963, "B-", 11563.897, 95.934133398, 3.599);
    static const Species Sr_96("Sr-96", "Sr", 20, 58, 38, 96, 8521.2041, "B-", 5411.738, 95.921719045, 9.089);
    static const Species Y_96("Y-96", "Y", 18, 57, 39, 96, 8569.4269, "B-", 7108.8741, 95.915909305, 6.521);
    static const Species Zr_96("Zr-96", "Zr", 16, 56, 40, 96, 8635.3283, "B-", 163.9704, 95.908277615, 0.122);
    static const Species Nb_96("Nb-96", "Nb", 14, 55, 41, 96, 8628.8868, "B-", 3192.059, 95.908101586, 0.157);
    static const Species Mo_96("Mo-96", "Mo", 12, 54, 42, 96, 8653.988, "B-", -2973.2411, 95.90467477, 0.128);
    static const Species Tc_96("Tc-96", "Tc", 10, 53, 43, 96, 8614.8673, "B-", 258.7369, 95.907866675, 5.524);
    static const Species Ru_96("Ru-96", "Ru", 8, 52, 44, 96, 8609.413, "B-", -6392.6529, 95.90758891, 0.182);
    static const Species Rh_96("Rh-96", "Rh", 6, 51, 45, 96, 8534.6735, "B-", -3504.3127, 95.914451705, 10.737);
    static const Species Pd_96("Pd-96", "Pd", 4, 50, 46, 96, 8490.0207, "B-", -11671.7741, 95.918213739, 4.502);
    static const Species Ag_96("Ag-96", "Ag", 2, 49, 47, 96, 8360.2903, "B-", -8940.0, 95.930743903, 96.708);
    static const Species Cd_96("Cd-96", "Cd", 0, 48, 48, 96, 8259.0, "B-", -17482.0, 95.940341, 440.0);
    static const Species In_96("In-96", "In", -2, 47, 49, 96, 8069.0, "B-", std::numeric_limits<double>::quiet_NaN(), 95.959109, 537.0);
    static const Species Br_97("Br-97", "Br", 27, 62, 35, 97, 8140.0, "B-", 13423.0, 96.963499, 429.0);
    static const Species Kr_97("Kr-97", "Kr", 25, 61, 36, 97, 8269.8645, "B-", 11095.646, 96.949088782, 140.0);
    static const Species Rb_97("Rb-97", "Rb", 23, 60, 37, 97, 8376.1872, "B-", 10061.5295, 96.937177117, 2.052);
    static const Species Sr_97("Sr-97", "Sr", 21, 59, 38, 97, 8471.8489, "B-", 7534.7807, 96.926375621, 3.633);
    static const Species Y_97("Y-97", "Y", 19, 58, 39, 97, 8541.4616, "B-", 6821.2382, 96.918286702, 7.201);
    static const Species Zr_97("Zr-97", "Zr", 17, 57, 40, 97, 8603.7182, "B-", 2666.1038, 96.910963802, 0.13);
    static const Species Nb_97("Nb-97", "Nb", 15, 56, 41, 97, 8623.1384, "B-", 1941.9038, 96.908101622, 4.556);
    static const Species Mo_97("Mo-97", "Mo", 13, 55, 42, 97, 8635.0926, "B-", -320.264, 96.906016903, 0.176);
    static const Species Tc_97("Tc-97", "Tc", 11, 54, 43, 97, 8623.7254, "B-", -1103.8722, 96.90636072, 4.42);
    static const Species Ru_97("Ru-97", "Ru", 9, 53, 44, 97, 8604.2799, "B-", -3523.0, 96.907545776, 2.965);
    static const Species Rh_97("Rh-97", "Rh", 7, 52, 45, 97, 8559.8949, "B-", -4791.7118, 96.911327872, 38.071);
    static const Species Pd_97("Pd-97", "Pd", 5, 51, 46, 97, 8502.4303, "B-", -6901.8255, 96.916471985, 5.2);
    static const Species Ag_97("Ag-97", "Ag", 3, 50, 47, 97, 8423.2121, "B-", -10170.0, 96.9238814, 12.9);
    static const Species Cd_97("Cd-97", "Cd", 1, 49, 48, 97, 8310.3013, "B-", -13344.0, 96.934799343, 451.073);
    static const Species In_97("In-97", "In", -1, 48, 49, 97, 8165.0, "B-", std::numeric_limits<double>::quiet_NaN(), 96.949125, 430.0);
    static const Species Br_98("Br-98", "Br", 28, 63, 35, 98, 8078.0, "B-", 16070.0, 97.969887, 429.0);
    static const Species Kr_98("Kr-98", "Kr", 26, 62, 36, 98, 8234.0, "B-", 10249.0, 97.952635, 322.0);
    static const Species Rb_98("Rb-98", "Rb", 24, 61, 37, 98, 8330.7294, "B-", 12053.2361, 97.941632317, 17.265);
    static const Species Sr_98("Sr-98", "Sr", 22, 60, 38, 98, 8445.7385, "B-", 5866.3591, 97.928692636, 3.463);
    static const Species Y_98("Y-98", "Y", 20, 59, 39, 98, 8497.6162, "B-", 8993.0098, 97.922394841, 8.501);
    static const Species Zr_98("Zr-98", "Zr", 18, 58, 40, 98, 8581.3984, "B-", 2242.8547, 97.912740448, 9.065);
    static const Species Nb_98("Nb-98", "Nb", 16, 57, 41, 98, 8596.3016, "B-", 4591.3681, 97.910332645, 5.369);
    static const Species Mo_98("Mo-98", "Mo", 14, 56, 42, 98, 8635.1691, "B-", -1683.7664, 97.905403609, 0.186);
    static const Species Tc_98("Tc-98", "Tc", 12, 55, 43, 98, 8610.0047, "B-", 1792.6575, 97.907211206, 3.628);
    static const Species Ru_98("Ru-98", "Ru", 10, 54, 44, 98, 8620.314, "B-", -5049.6529, 97.905286709, 6.937);
    static const Species Rh_98("Rh-98", "Rh", 8, 53, 45, 98, 8560.8038, "B-", -1854.2331, 97.910707734, 12.782);
    static const Species Pd_98("Pd-98", "Pd", 6, 52, 46, 98, 8533.8999, "B-", -8254.5607, 97.912698335, 5.09);
    static const Species Ag_98("Ag-98", "Ag", 4, 51, 47, 98, 8441.6866, "B-", -5430.0, 97.92155997, 35.327);
    static const Species Cd_98("Cd-98", "Cd", 2, 50, 48, 98, 8378.2953, "B-", -13730.0, 97.927389315, 55.605);
    static const Species In_98("In-98", "In", 0, 49, 49, 98, 8230.0, "B-", std::numeric_limits<double>::quiet_NaN(), 97.942129, 327.0);
    static const Species Kr_99("Kr-99", "Kr", 27, 63, 36, 99, 8175.0, "B-", 12721.0, 98.958776, 429.0);
    static const Species Rb_99("Rb-99", "Rb", 25, 62, 37, 99, 8295.301, "B-", 11397.3767, 98.94511919, 4.327);
    static const Species Sr_99("Sr-99", "Sr", 23, 61, 38, 99, 8402.5235, "B-", 8125.2037, 98.932883604, 5.085);
    static const Species Y_99("Y-99", "Y", 21, 60, 39, 99, 8476.6938, "B-", 6972.9398, 98.924160839, 7.101);
    static const Species Zr_99("Zr-99", "Zr", 19, 59, 40, 99, 8539.225, "B-", 4718.6736, 98.916675081, 11.271);
    static const Species Nb_99("Nb-99", "Nb", 17, 58, 41, 99, 8578.9859, "B-", 3634.7623, 98.911609377, 12.886);
    static const Species Mo_99("Mo-99", "Mo", 15, 57, 42, 99, 8607.7982, "B-", 1357.7631, 98.907707299, 0.245);
    static const Species Tc_99("Tc-99", "Tc", 13, 56, 43, 99, 8613.6105, "B-", 297.5156, 98.906249681, 0.974);
    static const Species Ru_99("Ru-99", "Ru", 11, 55, 44, 99, 8608.7132, "B-", -2040.8632, 98.905930284, 0.368);
    static const Species Rh_99("Rh-99", "Rh", 9, 54, 45, 99, 8580.1959, "B-", -3401.6603, 98.908121241, 20.881);
    static const Species Pd_99("Pd-99", "Pd", 7, 53, 46, 99, 8537.9332, "B-", -5470.3785, 98.911773073, 5.482);
    static const Species Ag_99("Ag-99", "Ag", 5, 52, 47, 99, 8474.7744, "B-", -6781.3511, 98.917645766, 6.725);
    static const Species Cd_99("Cd-99", "Cd", 3, 51, 48, 99, 8398.3734, "B-", -8555.0, 98.924925845, 1.7);
    static const Species In_99("In-99", "In", 1, 50, 49, 99, 8304.0, "B-", -13400.0, 98.93411, 320.0);
    static const Species Sn_99("Sn-99", "Sn", -1, 49, 50, 99, 8161.0, "B-", std::numeric_limits<double>::quiet_NaN(), 98.948495, 625.0);
    static const Species Kr_100("Kr-100", "Kr", 28, 64, 36, 100, 8134.0, "B-", 11796.0, 99.962995, 429.0);
    static const Species Rb_100("Rb-100", "Rb", 26, 63, 37, 100, 8244.5085, "B-", 13551.6204, 99.950331532, 14.089);
    static const Species Sr_100("Sr-100", "Sr", 24, 62, 38, 100, 8372.2012, "B-", 7503.7365, 99.93578327, 7.426);
    static const Species Y_100("Y-100", "Y", 22, 61, 39, 100, 8439.4151, "B-", 9051.4949, 99.927727678, 12.0);
    static const Species Zr_100("Zr-100", "Zr", 20, 60, 40, 100, 8522.1066, "B-", 3418.5098, 99.918010499, 8.742);
    static const Species Nb_100("Nb-100", "Nb", 18, 59, 41, 100, 8548.4683, "B-", 6401.7829, 99.914340578, 8.562);
    static const Species Mo_100("Mo-100", "Mo", 16, 58, 42, 100, 8604.6626, "B-", -172.0776, 99.907467982, 0.322);
    static const Species Tc_100("Tc-100", "Tc", 14, 57, 43, 100, 8595.1184, "B-", 3206.4401, 99.907652715, 1.45);
    static const Species Ru_100("Ru-100", "Ru", 12, 56, 44, 100, 8619.3593, "B-", -3636.2612, 99.90421046, 0.367);
    static const Species Rh_100("Rh-100", "Rh", 10, 55, 45, 100, 8575.1732, "B-", -378.4577, 99.908114147, 19.458);
    static const Species Pd_100("Pd-100", "Pd", 8, 54, 46, 100, 8563.5652, "B-", -7074.703, 99.908520438, 18.934);
    static const Species Ag_100("Ag-100", "Ag", 6, 53, 47, 100, 8484.9947, "B-", -3943.374, 99.916115443, 5.367);
    static const Species Cd_100("Cd-100", "Cd", 4, 52, 48, 100, 8437.7375, "B-", -10016.4492, 99.920348829, 1.8);
    static const Species In_100("In-100", "In", 2, 51, 49, 100, 8329.7495, "B-", -7030.0, 99.931101929, 2.4);
    static const Species Sn_100("Sn-100", "Sn", 0, 50, 50, 100, 8251.626, "B-", std::numeric_limits<double>::quiet_NaN(), 99.938648944, 257.661);
    static const Species Kr_101("Kr-101", "Kr", 29, 65, 36, 101, 8075.0, "B-", 13987.0, 100.969318, 537.0);
    static const Species Rb_101("Rb-101", "Rb", 27, 64, 37, 101, 8206.1753, "B-", 12757.4969, 100.954302, 22.0);
    static const Species Sr_101("Sr-101", "Sr", 25, 63, 38, 101, 8324.7411, "B-", 9729.8721, 100.940606264, 9.103);
    static const Species Y_101("Y-101", "Y", 23, 62, 39, 101, 8413.3305, "B-", 8106.2003, 100.930160817, 7.601);
    static const Species Zr_101("Zr-101", "Zr", 21, 61, 40, 101, 8485.8439, "B-", 5730.5011, 100.921458454, 8.944);
    static const Species Nb_101("Nb-101", "Nb", 19, 60, 41, 101, 8534.8355, "B-", 4628.4637, 100.915306508, 4.024);
    static const Species Mo_101("Mo-101", "Mo", 17, 59, 42, 101, 8572.9159, "B-", 2824.6411, 100.910337648, 0.331);
    static const Species Tc_101("Tc-101", "Tc", 15, 58, 43, 101, 8593.1366, "B-", 1613.52, 100.907305271, 25.768);
    static const Species Ru_101("Ru-101", "Ru", 13, 57, 44, 101, 8601.366, "B-", -545.6846, 100.905573086, 0.443);
    static const Species Rh_101("Rh-101", "Rh", 11, 56, 45, 101, 8588.2172, "B-", -1980.2833, 100.906158903, 6.27);
    static const Species Pd_101("Pd-101", "Pd", 9, 55, 46, 101, 8560.8644, "B-", -4097.7606, 100.908284824, 4.925);
    static const Species Ag_101("Ag-101", "Ag", 7, 54, 47, 101, 8512.5465, "B-", -5497.9186, 100.912683951, 5.193);
    static const Species Cd_101("Cd-101", "Cd", 5, 53, 48, 101, 8450.3657, "B-", -7291.5642, 100.918586209, 1.6);
    static const Species In_101("In-101", "In", 3, 52, 49, 101, 8370.426, "B-", -8239.277, 100.926414025, 12.519);
    static const Species Sn_101("Sn-101", "Sn", 1, 51, 50, 101, 8281.103, "B-", std::numeric_limits<double>::quiet_NaN(), 100.935259252, 322.068);
    static const Species Rb_102("Rb-102", "Rb", 28, 65, 37, 102, 8152.7443, "B-", 14906.9991, 101.960008, 89.0);
    static const Species Sr_102("Sr-102", "Sr", 26, 64, 38, 102, 8291.2213, "B-", 9013.3301, 101.944004679, 72.0);
    static const Species Y_102("Y-102", "Y", 24, 63, 39, 102, 8371.9172, "B-", 10408.7856, 101.934328471, 4.381);
    static const Species Zr_102("Zr-102", "Zr", 22, 62, 40, 102, 8466.294, "B-", 4716.838, 101.923154181, 9.401);
    static const Species Nb_102("Nb-102", "Nb", 20, 61, 41, 102, 8504.8675, "B-", 7262.6008, 101.918090447, 2.695);
    static const Species Mo_102("Mo-102", "Mo", 18, 60, 42, 102, 8568.3994, "B-", 1012.0557, 101.910293725, 8.916);
    static const Species Tc_102("Tc-102", "Tc", 16, 59, 43, 102, 8570.6514, "B-", 4533.5134, 101.909207239, 9.84);
    static const Species Ru_102("Ru-102", "Ru", 14, 58, 44, 102, 8607.4275, "B-", -2323.1187, 101.904340312, 0.446);
    static const Species Rh_102("Rh-102", "Rh", 12, 57, 45, 102, 8576.9818, "B-", 1119.647, 101.906834282, 6.88);
    static const Species Pd_102("Pd-102", "Pd", 10, 56, 46, 102, 8580.2887, "B-", -5656.2615, 101.905632292, 0.449);
    static const Species Ag_102("Ag-102", "Ag", 8, 55, 47, 102, 8517.165, "B-", -2587.0, 101.911704538, 8.771);
    static const Species Cd_102("Cd-102", "Cd", 6, 54, 48, 102, 8484.1322, "B-", -8964.8059, 101.914481797, 1.784);
    static const Species In_102("In-102", "In", 4, 53, 49, 102, 8388.5719, "B-", -5760.0, 101.924105911, 4.909);
    static const Species Sn_102("Sn-102", "Sn", 2, 52, 50, 102, 8324.4313, "B-", -13835.0, 101.930289525, 107.466);
    static const Species Sb_102("Sb-102", "Sb", 0, 51, 51, 102, 8181.0, "B-", std::numeric_limits<double>::quiet_NaN(), 101.945142, 429.0);
    static const Species Rb_103("Rb-103", "Rb", 29, 66, 37, 103, 8112.0, "B-", 14120.0, 102.964401, 429.0);
    static const Species Sr_103("Sr-103", "Sr", 27, 65, 38, 103, 8242.0, "B-", 11177.0, 102.949243, 215.0);
    static const Species Y_103("Y-103", "Y", 25, 64, 39, 103, 8342.6336, "B-", 9351.96, 102.937243796, 12.029);
    static const Species Zr_103("Zr-103", "Zr", 23, 63, 40, 103, 8425.8337, "B-", 7219.674, 102.927204054, 9.9);
    static const Species Nb_103("Nb-103", "Nb", 21, 62, 41, 103, 8488.332, "B-", 5925.6639, 102.919453416, 4.224);
    static const Species Mo_103("Mo-103", "Mo", 19, 61, 42, 103, 8538.2672, "B-", 3649.5889, 102.913091954, 9.9);
    static const Species Tc_103("Tc-103", "Tc", 17, 60, 43, 103, 8566.1045, "B-", 2663.2474, 102.90917396, 10.531);
    static const Species Ru_103("Ru-103", "Ru", 15, 59, 44, 103, 8584.3656, "B-", 764.5378, 102.906314846, 0.473);
    static const Species Rh_103("Rh-103", "Rh", 13, 58, 45, 103, 8584.1927, "B-", -574.7252, 102.905494081, 2.47);
    static const Species Pd_103("Pd-103", "Pd", 11, 57, 46, 103, 8571.0173, "B-", -2654.2778, 102.906111074, 0.942);
    static const Species Ag_103("Ag-103", "Ag", 9, 56, 47, 103, 8537.652, "B-", -4151.0761, 102.908960558, 4.4);
    static const Species Cd_103("Cd-103", "Cd", 7, 55, 48, 103, 8489.7547, "B-", -6019.2293, 102.913416922, 1.943);
    static const Species In_103("In-103", "In", 5, 54, 49, 103, 8423.7199, "B-", -7540.0, 102.91987883, 9.64);
    static const Species Sn_103("Sn-103", "Sn", 3, 53, 50, 103, 8343.0, "B-", -10422.0, 102.927973, 108.0);
    static const Species Sb_103("Sb-103", "Sb", 1, 52, 51, 103, 8234.0, "B-", std::numeric_limits<double>::quiet_NaN(), 102.939162, 322.0);
    static const Species Rb_104("Rb-104", "Rb", 30, 67, 37, 104, 8057.0, "B-", 16310.0, 103.970531, 537.0);
    static const Species Sr_104("Sr-104", "Sr", 28, 66, 38, 104, 8206.0, "B-", 10320.0, 103.953022, 322.0);
    static const Species Y_104("Y-104", "Y", 26, 65, 39, 104, 8298.0, "B-", 11638.0, 103.941943, 215.0);
    static const Species Zr_104("Zr-104", "Zr", 24, 64, 40, 104, 8402.3159, "B-", 6093.3367, 103.929449193, 10.0);
    static const Species Nb_104("Nb-104", "Nb", 22, 63, 41, 104, 8453.3832, "B-", 8532.7512, 103.922907728, 1.915);
    static const Species Mo_104("Mo-104", "Mo", 20, 62, 42, 104, 8527.9063, "B-", 2155.2212, 103.913747443, 9.566);
    static const Species Tc_104("Tc-104", "Tc", 18, 61, 43, 104, 8541.107, "B-", 5596.7945, 103.911433718, 26.716);
    static const Species Ru_104("Ru-104", "Ru", 16, 60, 44, 104, 8587.3998, "B-", -1136.4195, 103.905425312, 2.682);
    static const Species Rh_104("Rh-104", "Rh", 14, 59, 45, 104, 8568.9501, "B-", 2435.7789, 103.906645309, 2.471);
    static const Species Pd_104("Pd-104", "Pd", 12, 58, 46, 104, 8584.8485, "B-", -4278.6529, 103.904030393, 1.434);
    static const Species Ag_104("Ag-104", "Ag", 10, 57, 47, 104, 8536.185, "B-", -1148.0787, 103.908623715, 4.527);
    static const Species Cd_104("Cd-104", "Cd", 8, 56, 48, 104, 8517.6232, "B-", -7785.7166, 103.909856228, 1.795);
    static const Species In_104("In-104", "In", 6, 55, 49, 104, 8435.238, "B-", -4555.6174, 103.918214538, 6.2);
    static const Species Sn_104("Sn-104", "Sn", 4, 54, 50, 104, 8383.9114, "B-", -12332.0, 103.923105195, 6.167);
    static const Species Sb_104("Sb-104", "Sb", 2, 53, 51, 104, 8258.0, "B-", -9668.0, 103.936344, 109.0);
    static const Species Te_104("Te-104", "Te", 0, 52, 52, 104, 8157.3256, "B-", std::numeric_limits<double>::quiet_NaN(), 103.946723408, 340.967);
    static const Species Sr_105("Sr-105", "Sr", 29, 67, 38, 105, 8152.0, "B-", 12380.0, 104.959001, 537.0);
    static const Species Y_105("Y-105", "Y", 27, 66, 39, 105, 8262.0, "B-", 10888.0, 104.945711, 429.0);
    static const Species Zr_105("Zr-105", "Zr", 25, 65, 40, 105, 8358.598, "B-", 8457.2728, 104.934021832, 13.0);
    static const Species Nb_105("Nb-105", "Nb", 23, 64, 41, 105, 8431.6925, "B-", 7415.2411, 104.924942577, 4.324);
    static const Species Mo_105("Mo-105", "Mo", 21, 63, 42, 105, 8494.863, "B-", 4955.5157, 104.91698198899999, 9.721);
    static const Species Tc_105("Tc-105", "Tc", 19, 62, 43, 105, 8534.6074, "B-", 3648.2396, 104.911662024, 37.856);
    static const Species Ru_105("Ru-105", "Ru", 17, 61, 44, 105, 8561.9016, "B-", 1916.7271, 104.907745478, 2.683);
    static const Species Rh_105("Rh-105", "Rh", 15, 60, 45, 105, 8572.7053, "B-", 566.6347, 104.905687787, 2.685);
    static const Species Pd_105("Pd-105", "Pd", 13, 59, 46, 105, 8570.6509, "B-", -1347.0564, 104.905079479, 1.222);
    static const Species Ag_105("Ag-105", "Ag", 11, 58, 47, 105, 8550.3708, "B-", -2736.9989, 104.906525604, 4.877);
    static const Species Cd_105("Cd-105", "Cd", 9, 57, 48, 105, 8516.8532, "B-", -4693.2673, 104.909463893, 1.494);
    static const Species In_105("In-105", "In", 7, 56, 49, 105, 8464.7045, "B-", -6302.5807, 104.914502322, 11.0);
    static const Species Sn_105("Sn-105", "Sn", 5, 55, 50, 105, 8397.229, "B-", -9322.5103, 104.921268421, 4.263);
    static const Species Sb_105("Sb-105", "Sb", 3, 54, 51, 105, 8300.9923, "B-", -11203.9825, 104.931276547, 23.431);
    static const Species Te_105("Te-105", "Te", 1, 53, 52, 105, 8186.8368, "B-", std::numeric_limits<double>::quiet_NaN(), 104.943304516, 322.084);
    static const Species Sr_106("Sr-106", "Sr", 30, 68, 38, 106, 8114.0, "B-", 11490.0, 105.963177, 644.0);
    static const Species Y_106("Y-106", "Y", 28, 67, 39, 106, 8215.0, "B-", 12959.0, 105.950842, 537.0);
    static const Species Zr_106("Zr-106", "Zr", 26, 66, 40, 106, 8330.0, "B-", 7453.0, 105.93693, 215.0);
    static const Species Nb_106("Nb-106", "Nb", 24, 65, 41, 106, 8393.2657, "B-", 9925.3249, 105.928928505, 1.52);
    static const Species Mo_106("Mo-106", "Mo", 22, 64, 42, 106, 8479.5202, "B-", 3648.2494, 105.918273231, 9.801);
    static const Species Tc_106("Tc-106", "Tc", 20, 63, 43, 106, 8506.557, "B-", 6547.0, 105.914356674, 13.15);
    static const Species Ru_106("Ru-106", "Ru", 18, 62, 44, 106, 8560.9406, "B-", 39.4038, 105.907328181, 5.787);
    static const Species Rh_106("Rh-106", "Rh", 16, 61, 45, 106, 8553.9317, "B-", 3544.8865, 105.907285879, 5.786);
    static const Species Pd_106("Pd-106", "Pd", 14, 60, 46, 106, 8579.9934, "B-", -2965.1434, 105.903480287, 1.186);
    static const Species Ag_106("Ag-106", "Ag", 12, 59, 47, 106, 8544.6397, "B-", 189.7534, 105.906663499, 3.237);
    static const Species Cd_106("Cd-106", "Cd", 10, 58, 48, 106, 8539.0492, "B-", -6524.0031, 105.906459791, 1.184);
    static const Species In_106("In-106", "In", 8, 57, 49, 106, 8470.1213, "B-", -3254.4521, 105.913463596, 13.125);
    static const Species Sn_106("Sn-106", "Sn", 6, 56, 50, 106, 8432.0383, "B-", -10880.3964, 105.916957394, 5.465);
    static const Species Sb_106("Sb-106", "Sb", 4, 55, 51, 106, 8322.0124, "B-", -8253.5423, 105.928637979, 8.0);
    static const Species Te_106("Te-106", "Te", 2, 54, 52, 106, 8236.7682, "B-", -14920.0, 105.937498521, 107.934);
    static const Species I_106("I-106", "I", 0, 53, 53, 106, 8089.0, "B-", std::numeric_limits<double>::quiet_NaN(), 105.953516, 429.0);
    static const Species Sr_107("Sr-107", "Sr", 31, 69, 38, 107, 8057.0, "B-", 13720.0, 106.969672, 751.0);
    static const Species Y_107("Y-107", "Y", 29, 68, 39, 107, 8178.0, "B-", 12050.0, 106.954943, 537.0);
    static const Species Zr_107("Zr-107", "Zr", 27, 67, 40, 107, 8284.0, "B-", 9704.0, 106.942007, 322.0);
    static const Species Nb_107("Nb-107", "Nb", 25, 66, 41, 107, 8367.0898, "B-", 8821.1703, 106.931589685, 8.612);
    static const Species Mo_107("Mo-107", "Mo", 23, 65, 42, 107, 8442.219, "B-", 6204.9921, 106.92211977, 9.901);
    static const Species Tc_107("Tc-107", "Tc", 21, 64, 43, 107, 8492.8979, "B-", 5112.5985, 106.915458437, 9.31);
    static const Species Ru_107("Ru-107", "Ru", 19, 63, 44, 107, 8533.3676, "B-", 3001.1457, 106.909969837, 9.31);
    static const Species Rh_107("Rh-107", "Rh", 17, 62, 45, 107, 8554.104, "B-", 1508.9427, 106.906747975, 12.937);
    static const Species Pd_107("Pd-107", "Pd", 15, 61, 46, 107, 8560.8946, "B-", 34.0458, 106.905128058, 1.289);
    static const Species Ag_107("Ag-107", "Ag", 13, 60, 47, 107, 8553.9012, "B-", -1416.3741, 106.905091509, 2.556);
    static const Species Cd_107("Cd-107", "Cd", 11, 59, 48, 107, 8533.3524, "B-", -3423.6586, 106.906612049, 1.782);
    static const Species In_107("In-107", "In", 9, 58, 49, 107, 8494.0439, "B-", -5054.4281, 106.910287497, 10.363);
    static const Species Sn_107("Sn-107", "Sn", 7, 57, 50, 107, 8439.4946, "B-", -7858.9903, 106.915713649, 5.7);
    static const Species Sb_107("Sb-107", "Sb", 5, 56, 51, 107, 8358.7344, "B-", -9996.0, 106.924150621, 4.452);
    static const Species Te_107("Te-107", "Te", 3, 55, 52, 107, 8258.0, "B-", -11227.0, 106.934882, 108.0);
    static const Species I_107("I-107", "I", 1, 54, 53, 107, 8146.0, "B-", std::numeric_limits<double>::quiet_NaN(), 106.946935, 322.0);
    static const Species Y_108("Y-108", "Y", 30, 69, 39, 108, 8129.0, "B-", 14170.0, 107.960515, 644.0);
    static const Species Zr_108("Zr-108", "Zr", 28, 68, 40, 108, 8253.0, "B-", 8595.0, 107.945303, 429.0);
    static const Species Nb_108("Nb-108", "Nb", 26, 67, 41, 108, 8325.6604, "B-", 11204.0998, 107.936075604, 8.844);
    static const Species Mo_108("Mo-108", "Mo", 24, 66, 42, 108, 8422.1581, "B-", 5173.533, 107.924047508, 9.901);
    static const Species Tc_108("Tc-108", "Tc", 22, 65, 43, 108, 8462.8172, "B-", 7738.5736, 107.918493493, 9.413);
    static const Species Ru_108("Ru-108", "Ru", 20, 64, 44, 108, 8527.2267, "B-", 1369.7517, 107.910185793, 9.318);
    static const Species Rh_108("Rh-108", "Rh", 18, 63, 45, 108, 8532.6657, "B-", 4493.0596, 107.908715304, 15.026);
    static const Species Pd_108("Pd-108", "Pd", 16, 62, 46, 108, 8567.0241, "B-", -1917.4238, 107.903891806, 1.189);
    static const Species Ag_108("Ag-108", "Ag", 14, 61, 47, 108, 8542.0262, "B-", 1645.6311, 107.905950245, 2.563);
    static const Species Cd_108("Cd-108", "Cd", 12, 60, 48, 108, 8550.0196, "B-", -5132.5944, 107.904183588, 1.205);
    static const Species In_108("In-108", "In", 10, 59, 49, 108, 8495.2516, "B-", -2049.8794, 107.909693654, 9.276);
    static const Species Sn_108("Sn-108", "Sn", 8, 58, 50, 108, 8469.0273, "B-", -9624.6079, 107.91189429, 5.778);
    static const Species Sb_108("Sb-108", "Sb", 6, 57, 51, 108, 8372.6666, "B-", -6663.6646, 107.922226731, 5.9);
    static const Species Te_108("Te-108", "Te", 4, 56, 52, 108, 8303.7221, "B-", -13011.0, 107.929380469, 5.808);
    static const Species I_108("I-108", "I", 2, 55, 53, 108, 8176.0, "B-", -10139.0, 107.943348, 109.0);
    static const Species Xe_108("Xe-108", "Xe", 0, 54, 54, 108, 8074.8886, "B-", std::numeric_limits<double>::quiet_NaN(), 107.954232285, 407.406);
    static const Species Y_109("Y-109", "Y", 31, 70, 39, 109, 8089.0, "B-", 13250.0, 108.965131, 751.0);
    static const Species Zr_109("Zr-109", "Zr", 29, 69, 40, 109, 8204.0, "B-", 10960.0, 108.950907, 537.0);
    static const Species Nb_109("Nb-109", "Nb", 27, 68, 41, 109, 8297.1307, "B-", 9969.4851, 108.939141, 462.5);
    static const Species Mo_109("Mo-109", "Mo", 25, 67, 42, 109, 8381.4163, "B-", 7623.5438, 108.928438318, 12.0);
    static const Species Tc_109("Tc-109", "Tc", 23, 66, 43, 109, 8444.1796, "B-", 6455.6267, 108.920254107, 10.38);
    static const Species Ru_109("Ru-109", "Ru", 21, 65, 44, 109, 8496.228, "B-", 4260.7958, 108.913323707, 9.612);
    static const Species Rh_109("Rh-109", "Rh", 19, 64, 45, 109, 8528.1404, "B-", 2607.2327, 108.908749555, 4.336);
    static const Species Pd_109("Pd-109", "Pd", 17, 63, 46, 109, 8544.8825, "B-", 1112.9469, 108.905950576, 1.195);
    static const Species Ag_109("Ag-109", "Ag", 15, 62, 47, 109, 8547.9155, "B-", -215.1002, 108.904755778, 1.381);
    static const Species Cd_109("Cd-109", "Cd", 13, 61, 48, 109, 8538.7646, "B-", -2014.8047, 108.904986697, 1.649);
    static const Species In_109("In-109", "In", 11, 60, 49, 109, 8513.1027, "B-", -3859.3453, 108.907149679, 4.261);
    static const Species Sn_109("Sn-109", "Sn", 9, 59, 50, 109, 8470.5183, "B-", -6379.194, 108.911292857, 8.533);
    static const Species Sb_109("Sb-109", "Sb", 7, 58, 51, 109, 8404.8161, "B-", -8535.5871, 108.918141203, 5.652);
    static const Species Te_109("Te-109", "Te", 5, 57, 52, 109, 8319.3305, "B-", -10042.8941, 108.927304532, 4.704);
    static const Species I_109("I-109", "I", 3, 56, 53, 109, 8220.0164, "B-", -11502.9589, 108.938086022, 7.223);
    static const Species Xe_109("Xe-109", "Xe", 1, 55, 54, 109, 8107.3071, "B-", std::numeric_limits<double>::quiet_NaN(), 108.950434955, 322.178);
    static const Species Zr_110("Zr-110", "Zr", 30, 70, 40, 110, 8171.0, "B-", 10090.0, 109.954675, 537.0);
    static const Species Nb_110("Nb-110", "Nb", 28, 69, 41, 110, 8255.2607, "B-", 12225.9002, 109.943843, 900.0);
    static const Species Mo_110("Mo-110", "Mo", 26, 68, 42, 110, 8359.293, "B-", 6498.7491, 109.930717956, 26.0);
    static const Species Tc_110("Tc-110", "Tc", 24, 67, 43, 110, 8411.2603, "B-", 9038.0654, 109.923741263, 10.195);
    static const Species Ru_110("Ru-110", "Ru", 22, 66, 44, 110, 8486.3123, "B-", 2756.0638, 109.914038501, 9.58);
    static const Species Rh_110("Rh-110", "Rh", 20, 65, 45, 110, 8504.2551, "B-", 5502.2116, 109.911079745, 19.114);
    static const Species Pd_110("Pd-110", "Pd", 18, 64, 46, 110, 8547.163, "B-", -873.5982, 109.905172878, 0.657);
    static const Species Ag_110("Ag-110", "Ag", 16, 63, 47, 110, 8532.1089, "B-", 2890.6633, 109.906110724, 1.38);
    static const Species Cd_110("Cd-110", "Cd", 14, 62, 48, 110, 8551.2755, "B-", -3878.0, 109.90300747, 0.407);
    static const Species In_110("In-110", "In", 12, 61, 49, 110, 8508.9087, "B-", -627.9769, 109.907170674, 12.402);
    static const Species Sn_110("Sn-110", "Sn", 10, 60, 50, 110, 8496.0875, "B-", -8392.248, 109.907844835, 14.79);
    static const Species Sb_110("Sb-110", "Sb", 8, 59, 51, 110, 8412.6821, "B-", -5219.924, 109.916854283, 6.4);
    static const Species Te_110("Te-110", "Te", 6, 58, 52, 110, 8358.116, "B-", -11761.9766, 109.922458102, 7.058);
    static const Species I_110("I-110", "I", 4, 57, 53, 110, 8244.0767, "B-", -8545.2075, 109.935085102, 66.494);
    static const Species Xe_110("Xe-110", "Xe", 2, 56, 54, 110, 8159.2808, "B-", std::numeric_limits<double>::quiet_NaN(), 109.944258759, 108.415);
    static const Species Zr_111("Zr-111", "Zr", 31, 71, 40, 111, 8118.0, "B-", 12480.0, 110.960837, 644.0);
    static const Species Nb_111("Nb-111", "Nb", 29, 70, 41, 111, 8223.0, "B-", 10980.0, 110.947439, 322.0);
    static const Species Mo_111("Mo-111", "Mo", 27, 69, 42, 111, 8315.2932, "B-", 9084.862, 110.935651966, 13.503);
    static const Species Tc_111("Tc-111", "Tc", 25, 68, 43, 111, 8390.0906, "B-", 7760.65, 110.925898966, 11.359);
    static const Species Ru_111("Ru-111", "Ru", 23, 67, 44, 111, 8452.9582, "B-", 5518.5456, 110.917567566, 10.394);
    static const Species Rh_111("Rh-111", "Rh", 21, 66, 45, 111, 8495.6267, "B-", 3682.0153, 110.911643164, 7.356);
    static const Species Pd_111("Pd-111", "Pd", 19, 65, 46, 111, 8521.7498, "B-", 2229.5607, 110.907690358, 0.785);
    static const Species Ag_111("Ag-111", "Ag", 17, 64, 47, 111, 8534.7878, "B-", 1036.8, 110.905296827, 1.565);
    static const Species Cd_111("Cd-111", "Cd", 15, 63, 48, 111, 8537.0801, "B-", -860.1972, 110.904183776, 0.383);
    static const Species In_111("In-111", "In", 13, 62, 49, 111, 8522.2824, "B-", -2453.4692, 110.905107236, 3.675);
    static const Species Sn_111("Sn-111", "Sn", 11, 61, 50, 111, 8493.131, "B-", -5101.834, 110.907741143, 5.728);
    static const Species Sb_111("Sb-111", "Sb", 9, 60, 51, 111, 8440.1203, "B-", -7249.2597, 110.913218187, 9.5);
    static const Species Te_111("Te-111", "Te", 7, 59, 52, 111, 8367.7635, "B-", -8633.6922, 110.921000587, 6.9);
    static const Species I_111("I-111", "I", 5, 58, 53, 111, 8282.9343, "B-", -10434.0, 110.930269236, 5.103);
    static const Species Xe_111("Xe-111", "Xe", 3, 57, 54, 111, 8182.0, "B-", -11620.0, 110.94147, 124.0);
    static const Species Cs_111("Cs-111", "Cs", 1, 56, 55, 111, 8070.0, "B-", std::numeric_limits<double>::quiet_NaN(), 110.953945, 215.0);
    static const Species Zr_112("Zr-112", "Zr", 32, 72, 40, 112, 8081.0, "B-", 11650.0, 111.965196, 751.0);
    static const Species Nb_112("Nb-112", "Nb", 30, 71, 41, 112, 8178.0, "B-", 13410.0, 111.952689, 322.0);
    static const Species Mo_112("Mo-112", "Mo", 28, 70, 42, 112, 8291.0, "B-", 7779.0, 111.938293, 215.0);
    static const Species Tc_112("Tc-112", "Tc", 26, 69, 43, 112, 8353.6217, "B-", 10371.9409, 111.929941658, 5.92);
    static const Species Ru_112("Ru-112", "Ru", 24, 68, 44, 112, 8439.2431, "B-", 4100.179, 111.918806922, 10.305);
    static const Species Rh_112("Rh-112", "Rh", 22, 67, 45, 112, 8468.8666, "B-", 6589.9874, 111.914405199, 47.327);
    static const Species Pd_112("Pd-112", "Pd", 20, 66, 46, 112, 8520.7205, "B-", 262.6897, 111.907330557, 7.027);
    static const Species Ag_112("Ag-112", "Ag", 18, 65, 47, 112, 8516.0807, "B-", 3991.1283, 111.907048548, 2.6);
    static const Species Cd_112("Cd-112", "Cd", 16, 64, 48, 112, 8544.7306, "B-", -2584.7306, 111.902763896, 0.268);
    static const Species In_112("In-112", "In", 14, 63, 49, 112, 8514.6674, "B-", 664.9224, 111.905538718, 4.563);
    static const Species Sn_112("Sn-112", "Sn", 12, 62, 50, 112, 8513.6189, "B-", -7056.076, 111.904824894, 0.315);
    static const Species Sb_112("Sb-112", "Sb", 10, 61, 51, 112, 8443.633, "B-", -4031.455, 111.912399903, 19.14);
    static const Species Te_112("Te-112", "Te", 8, 60, 52, 112, 8400.6527, "B-", -10504.1795, 111.916727848, 9.0);
    static const Species I_112("I-112", "I", 6, 59, 53, 112, 8299.8801, "B-", -7036.991, 111.928004548, 11.0);
    static const Species Xe_112("Xe-112", "Xe", 4, 58, 54, 112, 8230.0646, "B-", -13612.0, 111.935559068, 8.891);
    static const Species Cs_112("Cs-112", "Cs", 2, 57, 55, 112, 8102.0, "B-", std::numeric_limits<double>::quiet_NaN(), 111.950172, 124.0);
    static const Species Zr_113("Zr-113", "Zr", 33, 73, 40, 113, 8027.0, "B-", 13870.0, 112.971723, 322.0);
    static const Species Nb_113("Nb-113", "Nb", 31, 72, 41, 113, 8143.0, "B-", 12440.0, 112.956833, 429.0);
    static const Species Mo_113("Mo-113", "Mo", 29, 71, 42, 113, 8246.0, "B-", 10162.0, 112.943478, 322.0);
    static const Species Tc_113("Tc-113", "Tc", 27, 70, 43, 113, 8329.4652, "B-", 9056.2674, 112.932569032, 3.6);
    static const Species Ru_113("Ru-113", "Ru", 25, 69, 44, 113, 8402.6857, "B-", 6899.1276, 112.922846729, 41.097);
    static const Species Rh_113("Rh-113", "Rh", 23, 68, 45, 113, 8456.8165, "B-", 4823.5559, 112.915440212, 7.656);
    static const Species Pd_113("Pd-113", "Pd", 21, 67, 46, 113, 8492.5795, "B-", 3436.3252, 112.910261912, 7.458);
    static const Species Ag_113("Ag-113", "Ag", 19, 66, 47, 113, 8516.066, "B-", 2016.4615, 112.906572865, 17.866);
    static const Species Cd_113("Cd-113", "Cd", 17, 65, 48, 113, 8526.9874, "B-", 323.837, 112.904408105, 0.262);
    static const Species In_113("In-113", "In", 15, 64, 49, 113, 8522.9297, "B-", -1038.9941, 112.904060451, 0.202);
    static const Species Sn_113("Sn-113", "Sn", 13, 63, 50, 113, 8506.8117, "B-", -3911.1637, 112.905175857, 1.69);
    static const Species Sb_113("Sb-113", "Sb", 11, 62, 51, 113, 8465.2762, "B-", -6069.9281, 112.909374664, 18.457);
    static const Species Te_113("Te-113", "Te", 9, 61, 52, 113, 8404.6366, "B-", -7227.521, 112.915891, 30.0);
    static const Species I_113("I-113", "I", 7, 60, 53, 113, 8333.7528, "B-", -8915.8902, 112.923650062, 8.6);
    static const Species Xe_113("Xe-113", "Xe", 5, 59, 54, 113, 8247.9277, "B-", -10439.0876, 112.933221663, 7.342);
    static const Species Cs_113("Cs-113", "Cs", 3, 58, 55, 113, 8148.623, "B-", -12055.0, 112.944428484, 9.207);
    static const Species Ba_113("Ba-113", "Ba", 1, 57, 56, 113, 8035.0, "B-", std::numeric_limits<double>::quiet_NaN(), 112.95737, 322.0);
    static const Species Nb_114("Nb-114", "Nb", 32, 73, 41, 114, 8097.0, "B-", 14720.0, 113.962469, 537.0);
    static const Species Mo_114("Mo-114", "Mo", 30, 72, 42, 114, 8219.0, "B-", 8920.0, 113.946666, 322.0);
    static const Species Tc_114("Tc-114", "Tc", 28, 71, 43, 114, 8290.2599, "B-", 11620.919, 113.93709, 465.0);
    static const Species Ru_114("Ru-114", "Ru", 26, 70, 44, 114, 8385.3351, "B-", 5489.0622, 113.92461443, 3.817);
    static const Species Rh_114("Rh-114", "Rh", 24, 69, 45, 114, 8426.6221, "B-", 7780.0712, 113.91872168, 76.824);
    static const Species Pd_114("Pd-114", "Pd", 22, 68, 46, 114, 8488.0056, "B-", 1440.4642, 113.91036943, 7.459);
    static const Species Ag_114("Ag-114", "Ag", 20, 67, 47, 114, 8493.7786, "B-", 5084.1233, 113.908823029, 4.9);
    static const Species Cd_114("Cd-114", "Cd", 18, 66, 48, 114, 8531.5135, "B-", -1445.1268, 113.903364998, 0.296);
    static const Species In_114("In-114", "In", 16, 65, 49, 114, 8511.9742, "B-", 1989.9281, 113.904916405, 0.323);
    static const Species Sn_114("Sn-114", "Sn", 14, 64, 50, 114, 8522.5671, "B-", -6063.1189, 113.90278013, 0.031);
    static const Species Sb_114("Sb-114", "Sb", 12, 63, 51, 114, 8462.5191, "B-", -2606.9398, 113.909289155, 21.226);
    static const Species Te_114("Te-114", "Te", 10, 62, 52, 114, 8432.7885, "B-", -9250.7417, 113.91208782, 26.224);
    static const Species I_114("I-114", "I", 8, 61, 53, 114, 8344.779, "B-", -5553.036, 113.9220189, 21.5);
    static const Species Xe_114("Xe-114", "Xe", 6, 60, 54, 114, 8289.2054, "B-", -12399.9706, 113.927980329, 12.0);
    static const Species Cs_114("Cs-114", "Cs", 4, 59, 55, 114, 8173.5711, "B-", -8780.4915, 113.941292244, 91.323);
    static const Species Ba_114("Ba-114", "Ba", 2, 58, 56, 114, 8089.6865, "B-", std::numeric_limits<double>::quiet_NaN(), 113.950718489, 110.227);
    static const Species Nb_115("Nb-115", "Nb", 33, 74, 41, 115, 8061.0, "B-", 13670.0, 114.966849, 537.0);
    static const Species Mo_115("Mo-115", "Mo", 31, 73, 42, 115, 8173.0, "B-", 11247.0, 114.952174, 429.0);
    static const Species Tc_115("Tc-115", "Tc", 29, 72, 43, 115, 8264.0, "B-", 10309.0, 114.9401, 210.0);
    static const Species Ru_115("Ru-115", "Ru", 27, 71, 44, 115, 8346.814, "B-", 8123.9327, 114.929033049, 27.016);
    static const Species Rh_115("Rh-115", "Rh", 25, 70, 45, 115, 8410.6538, "B-", 6196.5938, 114.920311649, 7.857);
    static const Species Pd_115("Pd-115", "Pd", 23, 69, 46, 115, 8457.7342, "B-", 4556.7647, 114.913659333, 14.543);
    static const Species Ag_115("Ag-115", "Ag", 21, 68, 47, 115, 8490.5553, "B-", 3101.893, 114.908767445, 19.611);
    static const Species Cd_115("Cd-115", "Cd", 19, 67, 48, 115, 8510.7252, "B-", 1451.8768, 114.905437426, 0.699);
    static const Species In_115("In-115", "In", 17, 66, 49, 115, 8516.5472, "B-", 497.4892, 114.903878772, 0.012);
    static const Species Sn_115("Sn-115", "Sn", 15, 65, 50, 115, 8514.0702, "B-", -3030.4336, 114.903344695, 0.016);
    static const Species Sb_115("Sb-115", "Sb", 13, 64, 51, 115, 8480.9156, "B-", -4940.6447, 114.906598, 17.203);
    static const Species Te_115("Te-115", "Te", 11, 63, 52, 115, 8431.1504, "B-", -5724.9628, 114.911902, 30.0);
    static const Species I_115("I-115", "I", 9, 62, 53, 115, 8374.5651, "B-", -7681.0475, 114.918048, 31.0);
    static const Species Xe_115("Xe-115", "Xe", 7, 61, 54, 115, 8300.9704, "B-", -8957.0, 114.926293943, 13.0);
    static const Species Cs_115("Cs-115", "Cs", 5, 60, 55, 115, 8216.0, "B-", -10779.0, 114.93591, 110.0);
    static const Species Ba_115("Ba-115", "Ba", 3, 59, 56, 115, 8116.0, "B-", std::numeric_limits<double>::quiet_NaN(), 114.947482, 215.0);
    static const Species Nb_116("Nb-116", "Nb", 34, 75, 41, 116, 8012.0, "B-", 15980.0, 115.972914, 322.0);
    static const Species Mo_116("Mo-116", "Mo", 32, 74, 42, 116, 8143.0, "B-", 10003.0, 115.955759, 537.0);
    static const Species Tc_116("Tc-116", "Tc", 30, 73, 43, 116, 8223.0, "B-", 12855.0, 115.94502, 320.0);
    static const Species Ru_116("Ru-116", "Ru", 28, 72, 44, 116, 8326.884, "B-", 6666.8252, 115.931219191, 4.0);
    static const Species Rh_116("Rh-116", "Rh", 26, 71, 45, 116, 8377.6123, "B-", 9095.2839, 115.92406206, 79.261);
    static const Species Pd_116("Pd-116", "Pd", 24, 70, 46, 116, 8449.2755, "B-", 2711.6378, 115.914297872, 7.659);
    static const Species Ag_116("Ag-116", "Ag", 22, 69, 47, 116, 8465.9073, "B-", 6169.8248, 115.911386809, 3.5);
    static const Species Cd_116("Cd-116", "Cd", 20, 68, 48, 116, 8512.3511, "B-", -462.7305, 115.90476323, 0.172);
    static const Species In_116("In-116", "In", 18, 67, 49, 116, 8501.6177, "B-", 3276.2204, 115.905259992, 0.236);
    static const Species Sn_116("Sn-116", "Sn", 16, 66, 50, 116, 8523.1166, "B-", -4703.9591, 115.901742825, 0.103);
    static const Species Sb_116("Sb-116", "Sb", 14, 65, 51, 116, 8475.8208, "B-", -1558.2272, 115.906792732, 5.533);
    static const Species Te_116("Te-116", "Te", 12, 64, 52, 116, 8455.6435, "B-", -7843.1388, 115.908465558, 25.986);
    static const Species I_116("I-116", "I", 10, 63, 53, 116, 8381.2858, "B-", -4373.7764, 115.916885513, 80.555);
    static const Species Xe_116("Xe-116", "Xe", 8, 62, 54, 116, 8336.8365, "B-", -11004.0, 115.921580955, 13.974);
    static const Species Cs_116("Cs-116", "Cs", 6, 61, 55, 116, 8235.0, "B-", -7663.0, 115.933395, 108.0);
    static const Species Ba_116("Ba-116", "Ba", 4, 60, 56, 116, 8162.0, "B-", -14330.0, 115.941621, 215.0);
    static const Species La_116("La-116", "La", 2, 59, 57, 116, 8032.0, "B-", std::numeric_limits<double>::quiet_NaN(), 115.957005, 345.0);
    static const Species Mo_117("Mo-117", "Mo", 33, 75, 42, 117, 8096.0, "B-", 12450.0, 116.961686, 537.0);
    static const Species Tc_117("Tc-117", "Tc", 31, 74, 43, 117, 8195.0, "B-", 11350.0, 116.94832, 429.0);
    static const Species Ru_117("Ru-117", "Ru", 29, 73, 44, 117, 8285.5625, "B-", 9406.8875, 116.936135, 465.0);
    static const Species Rh_117("Rh-117", "Rh", 27, 72, 45, 117, 8359.2766, "B-", 7527.1313, 116.926036291, 9.548);
    static const Species Pd_117("Pd-117", "Pd", 25, 71, 46, 117, 8416.9243, "B-", 5758.0284, 116.917955584, 7.788);
    static const Species Ag_117("Ag-117", "Ag", 23, 70, 47, 117, 8459.4515, "B-", 4236.479, 116.911774086, 14.57);
    static const Species Cd_117("Cd-117", "Cd", 21, 69, 48, 117, 8488.974, "B-", 2524.6381, 116.907226039, 1.087);
    static const Species In_117("In-117", "In", 19, 68, 49, 117, 8503.8653, "B-", 1454.7073, 116.904515729, 5.239);
    static const Species Sn_117("Sn-117", "Sn", 17, 67, 50, 117, 8509.612, "B-", -1758.1788, 116.902954036, 0.518);
    static const Species Sb_117("Sb-117", "Sb", 15, 66, 51, 117, 8487.8981, "B-", -3544.0634, 116.904841519, 9.057);
    static const Species Te_117("Te-117", "Te", 13, 65, 52, 117, 8450.9203, "B-", -4656.9321, 116.908646227, 14.444);
    static const Species I_117("I-117", "I", 11, 64, 53, 117, 8404.4307, "B-", -6253.2213, 116.913645649, 27.437);
    static const Species Xe_117("Xe-117", "Xe", 9, 63, 54, 117, 8344.2976, "B-", -7692.2462, 116.920358758, 11.141);
    static const Species Cs_117("Cs-117", "Cs", 7, 62, 55, 117, 8271.8652, "B-", -9035.1943, 116.928616723, 67.0);
    static const Species Ba_117("Ba-117", "Ba", 5, 61, 56, 117, 8187.9546, "B-", -11187.0, 116.938316403, 268.749);
    static const Species La_117("La-117", "La", 3, 60, 57, 117, 8086.0, "B-", std::numeric_limits<double>::quiet_NaN(), 116.950326, 215.0);
    static const Species Mo_118("Mo-118", "Mo", 34, 76, 42, 118, 8067.0, "B-", 10920.0, 117.965249, 537.0);
    static const Species Tc_118("Tc-118", "Tc", 32, 75, 43, 118, 8153.0, "B-", 13710.0, 117.953526, 429.0);
    static const Species Ru_118("Ru-118", "Ru", 30, 74, 44, 118, 8263.0, "B-", 7887.0, 117.938808, 215.0);
    static const Species Rh_118("Rh-118", "Rh", 28, 73, 45, 118, 8322.8539, "B-", 10501.5182, 117.930341116, 26.018);
    static const Species Pd_118("Pd-118", "Pd", 26, 72, 46, 118, 8405.2197, "B-", 4165.4444, 117.919067273, 2.677);
    static const Species Ag_118("Ag-118", "Ag", 24, 71, 47, 118, 8433.89, "B-", 7147.8469, 117.914595484, 2.7);
    static const Species Cd_118("Cd-118", "Cd", 22, 70, 48, 118, 8487.835, "B-", 526.5277, 117.906921956, 21.471);
    static const Species In_118("In-118", "In", 20, 69, 49, 118, 8485.667, "B-", 4424.6664, 117.906356705, 8.322);
    static const Species Sn_118("Sn-118", "Sn", 18, 68, 50, 118, 8516.5341, "B-", -3656.6393, 117.90160663, 0.536);
    static const Species Sb_118("Sb-118", "Sb", 16, 67, 51, 118, 8478.9156, "B-", -305.4459, 117.905532194, 3.237);
    static const Species Te_118("Te-118", "Te", 14, 66, 52, 118, 8469.697, "B-", -6719.7015, 117.905860104, 19.652);
    static const Species I_118("I-118", "I", 12, 65, 53, 118, 8406.1203, "B-", -2891.9893, 117.913074, 21.213);
    static const Species Xe_118("Xe-118", "Xe", 10, 64, 54, 118, 8374.9819, "B-", -9669.6905, 117.916178678, 11.141);
    static const Species Cs_118("Cs-118", "Cs", 8, 63, 55, 118, 8286.4053, "B-", -6210.0, 117.926559517, 13.69);
    static const Species Ba_118("Ba-118", "Ba", 6, 62, 56, 118, 8227.0, "B-", -12580.0, 117.933226, 215.0);
    static const Species La_118("La-118", "La", 4, 61, 57, 118, 8114.0, "B-", std::numeric_limits<double>::quiet_NaN(), 117.946731, 322.0);
    static const Species Mo_119("Mo-119", "Mo", 35, 77, 42, 119, 8019.0, "B-", 13590.0, 118.971465, 322.0);
    static const Species Tc_119("Tc-119", "Tc", 33, 76, 43, 119, 8126.0, "B-", 11910.0, 118.956876, 537.0);
    static const Species Ru_119("Ru-119", "Ru", 31, 75, 44, 119, 8220.0, "B-", 10743.0, 118.94409, 322.0);
    static const Species Rh_119("Rh-119", "Rh", 29, 74, 45, 119, 8303.3953, "B-", 8584.4751, 118.932556951, 10.0);
    static const Species Pd_119("Pd-119", "Pd", 27, 73, 46, 119, 8368.9594, "B-", 7238.4816, 118.923341138, 8.854);
    static const Species Ag_119("Ag-119", "Ag", 25, 72, 47, 119, 8423.2126, "B-", 5331.1799, 118.915570309, 15.783);
    static const Species Cd_119("Cd-119", "Cd", 23, 71, 48, 119, 8461.4381, "B-", 3721.7197, 118.909847052, 40.467);
    static const Species In_119("In-119", "In", 21, 70, 49, 119, 8486.1387, "B-", 2366.3263, 118.905851622, 7.847);
    static const Species Sn_119("Sn-119", "Sn", 19, 69, 50, 119, 8499.4494, "B-", -589.4452, 118.903311266, 0.778);
    static const Species Sb_119("Sb-119", "Sb", 17, 68, 51, 119, 8487.9218, "B-", -2293.0, 118.903944062, 7.512);
    static const Species Te_119("Te-119", "Te", 15, 67, 52, 119, 8462.0785, "B-", -3404.808, 118.906405699, 7.813);
    static const Species I_119("I-119", "I", 13, 66, 53, 119, 8426.8924, "B-", -4983.2433, 118.91006091, 23.302);
    static const Species Xe_119("Xe-119", "Xe", 11, 65, 54, 119, 8378.442, "B-", -6489.4269, 118.915410641, 11.141);
    static const Species Cs_119("Cs-119", "Cs", 9, 64, 55, 119, 8317.3347, "B-", -7714.9651, 118.922377327, 14.965);
    static const Species Ba_119("Ba-119", "Ba", 7, 63, 56, 119, 8245.9287, "B-", -9570.0, 118.930659683, 214.997);
    static const Species La_119("La-119", "La", 5, 62, 57, 119, 8159.0, "B-", -11199.0, 118.940934, 322.0);
    static const Species Ce_119("Ce-119", "Ce", 3, 61, 58, 119, 8058.0, "B-", std::numeric_limits<double>::quiet_NaN(), 118.952957, 537.0);
    static const Species Tc_120("Tc-120", "Tc", 34, 77, 43, 120, 8083.0, "B-", 14720.0, 119.962426, 537.0);
    static const Species Ru_120("Ru-120", "Ru", 32, 76, 44, 120, 8199.0, "B-", 8899.0, 119.946623, 429.0);
    static const Species Rh_120("Rh-120", "Rh", 30, 75, 45, 120, 8266.0, "B-", 11660.0, 119.937069, 215.0);
    static const Species Pd_120("Pd-120", "Pd", 28, 74, 46, 120, 8357.0817, "B-", 5371.9076, 119.924551745, 2.464);
    static const Species Ag_120("Ag-120", "Ag", 26, 73, 47, 120, 8395.3281, "B-", 8305.8535, 119.918784765, 4.8);
    static const Species Cd_120("Cd-120", "Cd", 24, 72, 48, 120, 8458.024, "B-", 1770.3754, 119.909868065, 4.0);
    static const Species In_120("In-120", "In", 22, 71, 49, 120, 8466.2575, "B-", 5370.0, 119.907967489, 42.953);
    static const Species Sn_120("Sn-120", "Sn", 20, 70, 50, 120, 8504.488, "B-", -2680.6076, 119.902202557, 0.987);
    static const Species Sb_120("Sb-120", "Sb", 18, 69, 51, 120, 8475.63, "B-", 945.0271, 119.905080308, 7.728);
    static const Species Te_120("Te-120", "Te", 16, 68, 52, 120, 8476.9857, "B-", -5615.0, 119.904065779, 1.88);
    static const Species I_120("I-120", "I", 14, 67, 53, 120, 8423.6745, "B-", -1574.726, 119.910093729, 16.212);
    static const Species Xe_120("Xe-120", "Xe", 12, 66, 54, 120, 8404.0322, "B-", -8283.7857, 119.911784267, 12.686);
    static const Species Cs_120("Cs-120", "Cs", 10, 65, 55, 120, 8328.4811, "B-", -5000.0, 119.920677277, 10.702);
    static const Species Ba_120("Ba-120", "Ba", 8, 64, 56, 120, 8280.2949, "B-", -11319.0, 119.926044997, 322.241);
    static const Species La_120("La-120", "La", 6, 63, 57, 120, 8179.0, "B-", -7840.0, 119.938196, 322.0);
    static const Species Ce_120("Ce-120", "Ce", 4, 62, 58, 120, 8108.0, "B-", std::numeric_limits<double>::quiet_NaN(), 119.946613, 537.0);
    static const Species Tc_121("Tc-121", "Tc", 35, 78, 43, 121, 8054.0, "B-", 13080.0, 120.96614, 537.0);
    static const Species Ru_121("Ru-121", "Ru", 33, 77, 44, 121, 8156.0, "B-", 11630.0, 120.952098, 429.0);
    static const Species Rh_121("Rh-121", "Rh", 31, 76, 45, 121, 8245.2397, "B-", 9932.203, 120.939613, 665.0);
    static const Species Pd_121("Pd-121", "Pd", 29, 75, 46, 121, 8320.8584, "B-", 8220.4934, 120.928950342, 3.6);
    static const Species Ag_121("Ag-121", "Ag", 27, 74, 47, 121, 8382.3306, "B-", 6671.0057, 120.920125279, 13.0);
    static const Species Cd_121("Cd-121", "Cd", 25, 73, 48, 121, 8430.9972, "B-", 4760.7564, 120.91296366, 2.085);
    static const Species In_121("In-121", "In", 23, 72, 49, 121, 8463.8767, "B-", 3362.0331, 120.907852778, 29.435);
    static const Species Sn_121("Sn-121", "Sn", 21, 71, 50, 121, 8485.1964, "B-", 402.5306, 120.904243488, 1.05);
    static const Species Sb_121("Sb-121", "Sb", 19, 70, 51, 121, 8482.0574, "B-", -1056.0462, 120.903811353, 2.69);
    static const Species Te_121("Te-121", "Te", 17, 69, 52, 121, 8466.8641, "B-", -2297.4615, 120.904945065, 27.734);
    static const Species I_121("I-121", "I", 15, 68, 53, 121, 8441.4111, "B-", -3764.6525, 120.907411492, 5.07);
    static const Species Xe_121("Xe-121", "Xe", 13, 67, 54, 121, 8403.8326, "B-", -5378.6549, 120.911453012, 10.995);
    static const Species Cs_121("Cs-121", "Cs", 11, 66, 55, 121, 8352.9152, "B-", -6357.4948, 120.917227235, 15.34);
    static const Species Ba_121("Ba-121", "Ba", 9, 65, 56, 121, 8293.9083, "B-", -8555.0, 120.924052286, 152.333);
    static const Species La_121("La-121", "La", 7, 64, 57, 121, 8217.0, "B-", -9500.0, 120.933236, 322.0);
    static const Species Ce_121("Ce-121", "Ce", 5, 63, 58, 121, 8132.0, "B-", -11139.0, 120.943435, 430.0);
    static const Species Pr_121("Pr-121", "Pr", 3, 62, 59, 121, 8033.0, "B-", std::numeric_limits<double>::quiet_NaN(), 120.955393, 537.0);
    static const Species Tc_122("Tc-122", "Tc", 36, 79, 43, 122, 8011.0, "B-", 15475.0, 121.97176, 322.0);
    static const Species Ru_122("Ru-122", "Ru", 34, 78, 44, 122, 8132.0, "B-", 10099.0, 121.955147, 537.0);
    static const Species Rh_122("Rh-122", "Rh", 32, 77, 45, 122, 8208.0, "B-", 12737.0, 121.944305, 322.0);
    static const Species Pd_122("Pd-122", "Pd", 30, 76, 46, 122, 8305.9755, "B-", 6489.9492, 121.930631693, 21.0);
    static const Species Ag_122("Ag-122", "Ag", 28, 75, 47, 122, 8352.7591, "B-", 9506.2662, 121.923664446, 41.0);
    static const Species Cd_122("Cd-122", "Cd", 26, 74, 48, 122, 8424.2667, "B-", 2958.9765, 121.91345905, 2.468);
    static const Species In_122("In-122", "In", 24, 73, 49, 122, 8442.1079, "B-", 6368.5921, 121.910282458, 53.741);
    static const Species Sn_122("Sn-122", "Sn", 22, 72, 50, 122, 8487.8968, "B-", -1605.7483, 121.903445494, 2.627);
    static const Species Sb_122("Sb-122", "Sb", 20, 71, 51, 122, 8468.3222, "B-", 1979.0772, 121.905169335, 2.687);
    static const Species Te_122("Te-122", "Te", 18, 70, 52, 122, 8478.1315, "B-", -4234.0, 121.903044708, 1.456);
    static const Species I_122("I-122", "I", 16, 69, 53, 122, 8437.0139, "B-", -724.2937, 121.907590094, 5.561);
    static const Species Xe_122("Xe-122", "Xe", 14, 68, 54, 122, 8424.6644, "B-", -7210.2195, 121.908367655, 11.928);
    static const Species Cs_122("Cs-122", "Cs", 12, 67, 55, 122, 8359.1515, "B-", -3535.817, 121.916108144, 36.164);
    static const Species Ba_122("Ba-122", "Ba", 10, 66, 56, 122, 8323.7567, "B-", -10066.0, 121.919904, 30.0);
    static const Species La_122("La-122", "La", 8, 65, 57, 122, 8235.0, "B-", -6669.0, 121.93071, 320.0);
    static const Species Ce_122("Ce-122", "Ce", 6, 64, 58, 122, 8174.0, "B-", -13094.0, 121.93787, 430.0);
    static const Species Pr_122("Pr-122", "Pr", 4, 63, 59, 122, 8060.0, "B-", std::numeric_limits<double>::quiet_NaN(), 121.951927, 537.0);
    static const Species Ru_123("Ru-123", "Ru", 35, 79, 44, 123, 8089.0, "B-", 12640.0, 122.960762, 537.0);
    static const Species Rh_123("Rh-123", "Rh", 33, 78, 45, 123, 8185.0, "B-", 11239.0, 122.947192, 429.0);
    static const Species Pd_123("Pd-123", "Pd", 31, 77, 46, 123, 8270.0318, "B-", 9138.8323, 122.935126, 847.5);
    static const Species Ag_123("Ag-123", "Ag", 29, 76, 47, 123, 8337.9707, "B-", 7845.6026, 122.92531506, 35.0);
    static const Species Cd_123("Cd-123", "Cd", 27, 75, 48, 123, 8395.3955, "B-", 6014.8503, 122.91689246, 2.894);
    static const Species In_123("In-123", "In", 25, 74, 49, 123, 8437.9362, "B-", 4385.6489, 122.910435252, 21.29);
    static const Species Sn_123("Sn-123", "Sn", 23, 73, 50, 123, 8467.2313, "B-", 1408.2079, 122.905727065, 2.661);
    static const Species Sb_123("Sb-123", "Sb", 21, 72, 51, 123, 8472.3196, "B-", -51.9128, 122.904215292, 1.456);
    static const Species Te_123("Te-123", "Te", 19, 71, 52, 123, 8465.537, "B-", -1228.3898, 122.904271022, 1.454);
    static const Species I_123("I-123", "I", 17, 70, 53, 123, 8449.1896, "B-", -2694.3302, 122.905589753, 3.956);
    static const Species Xe_123("Xe-123", "Xe", 15, 69, 54, 123, 8420.9239, "B-", -4204.6012, 122.908482235, 10.234);
    static const Species Cs_123("Cs-123", "Cs", 13, 68, 55, 123, 8380.3796, "B-", -5388.6934, 122.91299606, 13.0);
    static const Species Ba_123("Ba-123", "Ba", 11, 67, 56, 123, 8330.2086, "B-", -7004.0, 122.91878106, 13.0);
    static const Species La_123("La-123", "La", 9, 66, 57, 123, 8267.0, "B-", -8365.0, 122.9263, 210.0);
    static const Species Ce_123("Ce-123", "Ce", 7, 65, 58, 123, 8193.0, "B-", -10056.0, 122.93528, 320.0);
    static const Species Pr_123("Pr-123", "Pr", 5, 64, 59, 123, 8104.0, "B-", std::numeric_limits<double>::quiet_NaN(), 122.946076, 429.0);
    static const Species Ru_124("Ru-124", "Ru", 36, 80, 44, 124, 8065.0, "B-", 11120.0, 123.96394, 644.0);
    static const Species Rh_124("Rh-124", "Rh", 34, 79, 45, 124, 8148.0, "B-", 13690.0, 123.952002, 429.0);
    static const Species Pd_124("Pd-124", "Pd", 32, 78, 46, 124, 8252.0, "B-", 7830.0, 123.937305, 322.0);
    static const Species Ag_124("Ag-124", "Ag", 30, 77, 47, 124, 8308.8958, "B-", 10469.4858, 123.928899227, 270.0);
    static const Species Cd_124("Cd-124", "Cd", 28, 76, 48, 124, 8387.0179, "B-", 4168.342, 123.917659772, 2.8);
    static const Species In_124("In-124", "In", 26, 75, 49, 124, 8414.3243, "B-", 7363.697, 123.913184873, 32.808);
    static const Species Sn_124("Sn-124", "Sn", 24, 74, 50, 124, 8467.3997, "B-", -612.4067, 123.905279619, 1.41);
    static const Species Sb_124("Sb-124", "Sb", 22, 73, 51, 124, 8456.1517, "B-", 2905.073, 123.905937065, 1.457);
    static const Species Te_124("Te-124", "Te", 20, 72, 52, 124, 8473.2705, "B-", -3159.587, 123.902818341, 1.451);
    static const Species I_124("I-124", "I", 18, 71, 53, 124, 8441.4807, "B-", 302.8501, 123.906210297, 2.467);
    static const Species Xe_124("Xe-124", "Xe", 16, 70, 54, 124, 8437.6138, "B-", -5926.3445, 123.905885174, 1.457);
    static const Species Cs_124("Cs-124", "Cs", 14, 69, 55, 124, 8383.5114, "B-", -2651.2748, 123.912247366, 9.823);
    static const Species Ba_124("Ba-124", "Ba", 12, 68, 56, 124, 8355.8209, "B-", -8831.1685, 123.915093627, 13.416);
    static const Species La_124("La-124", "La", 10, 67, 57, 124, 8278.2926, "B-", -5343.0, 123.924574275, 60.836);
    static const Species Ce_124("Ce-124", "Ce", 8, 66, 58, 124, 8229.0, "B-", -11765.0, 123.93031, 320.0);
    static const Species Pr_124("Pr-124", "Pr", 6, 65, 59, 124, 8128.0, "B-", -8321.0, 123.94294, 430.0);
    static const Species Nd_124("Nd-124", "Nd", 4, 64, 60, 124, 8054.0, "B-", std::numeric_limits<double>::quiet_NaN(), 123.951873, 537.0);
    static const Species Ru_125("Ru-125", "Ru", 37, 81, 44, 125, 8023.0, "B-", 13460.0, 124.969544, 322.0);
    static const Species Rh_125("Rh-125", "Rh", 35, 80, 45, 125, 8124.0, "B-", 12130.0, 124.955094, 537.0);
    static const Species Pd_125("Pd-125", "Pd", 33, 79, 46, 125, 8215.0, "B-", 10560.0, 124.942072, 429.0);
    static const Species Ag_125("Ag-125", "Ag", 31, 78, 47, 125, 8293.3151, "B-", 8828.1511, 124.930735, 465.0);
    static const Species Cd_125("Cd-125", "Cd", 29, 77, 48, 125, 8357.6815, "B-", 7064.2177, 124.92125759, 3.1);
    static const Species In_125("In-125", "In", 27, 76, 49, 125, 8407.9365, "B-", 5481.3495, 124.913673841, 1.9);
    static const Species Sn_125("Sn-125", "Sn", 25, 75, 50, 125, 8445.5285, "B-", 2361.4366, 124.90778937, 1.426);
    static const Species Sb_125("Sb-125", "Sb", 23, 74, 51, 125, 8458.1612, "B-", 766.7, 124.905254264, 2.7);
    static const Species Te_125("Te-125", "Te", 21, 73, 52, 125, 8458.0361, "B-", -185.77, 124.904431178, 1.451);
    static const Species I_125("I-125", "I", 19, 72, 53, 125, 8450.2911, "B-", -1636.6632, 124.90463061, 1.452);
    static const Species Xe_125("Xe-125", "Xe", 17, 71, 54, 125, 8430.939, "B-", -3109.6184, 124.90638764, 1.518);
    static const Species Cs_125("Cs-125", "Cs", 15, 70, 55, 125, 8399.8033, "B-", -4420.7663, 124.909725953, 8.304);
    static const Species Ba_125("Ba-125", "Ba", 13, 69, 56, 125, 8358.1784, "B-", -5909.4836, 124.91447184, 11.8);
    static const Species La_125("La-125", "La", 11, 68, 57, 125, 8304.6438, "B-", -7102.0, 124.920815931, 27.909);
    static const Species Ce_125("Ce-125", "Ce", 9, 67, 58, 125, 8242.0, "B-", -8587.0, 124.92844, 210.0);
    static const Species Pr_125("Pr-125", "Pr", 7, 66, 59, 125, 8167.0, "B-", -10001.0, 124.937659, 322.0);
    static const Species Nd_125("Nd-125", "Nd", 5, 65, 60, 125, 8080.0, "B-", std::numeric_limits<double>::quiet_NaN(), 124.948395, 429.0);
    static const Species Rh_126("Rh-126", "Rh", 36, 81, 45, 126, 8087.0, "B-", 14590.0, 125.960064, 537.0);
    static const Species Pd_126("Pd-126", "Pd", 34, 80, 46, 126, 8197.0, "B-", 8930.0, 125.944401, 429.0);
    static const Species Ag_126("Ag-126", "Ag", 32, 79, 47, 126, 8261.0, "B-", 11535.0, 125.934814, 215.0);
    static const Species Cd_126("Cd-126", "Cd", 30, 78, 48, 126, 8346.7393, "B-", 5553.65, 125.92243029, 2.473);
    static const Species In_126("In-126", "In", 28, 77, 49, 126, 8384.6067, "B-", 8205.7585, 125.916468202, 4.5);
    static const Species Sn_126("Sn-126", "Sn", 26, 76, 50, 126, 8443.5227, "B-", 378.0, 125.907658958, 11.473);
    static const Species Sb_126("Sb-126", "Sb", 24, 75, 51, 126, 8440.3136, "B-", 3671.0321, 125.907253158, 34.189);
    static const Species Te_126("Te-126", "Te", 22, 74, 52, 126, 8463.2397, "B-", -2153.6712, 125.903312144, 1.453);
    static const Species I_126("I-126", "I", 20, 73, 53, 126, 8439.9379, "B-", 1235.8904, 125.905624205, 4.055);
    static const Species Xe_126("Xe-126", "Xe", 18, 72, 54, 126, 8443.5375, "B-", -4795.7039, 125.904297422, 0.006);
    static const Species Cs_126("Cs-126", "Cs", 16, 71, 55, 126, 8399.2673, "B-", -1680.7697, 125.909445821, 11.12);
    static const Species Ba_126("Ba-126", "Ba", 14, 70, 56, 126, 8379.7187, "B-", -7696.4376, 125.911250202, 13.416);
    static const Species La_126("La-126", "La", 12, 69, 57, 126, 8312.4268, "B-", -4152.9106, 125.919512667, 97.163);
    static const Species Ce_126("Ce-126", "Ce", 10, 68, 58, 126, 8273.2581, "B-", -10497.0, 125.923971, 30.0);
    static const Species Pr_126("Pr-126", "Pr", 8, 67, 59, 126, 8184.0, "B-", -6943.0, 125.93524, 210.0);
    static const Species Nd_126("Nd-126", "Nd", 6, 66, 60, 126, 8122.0, "B-", -13631.0, 125.942694, 322.0);
    static const Species Pm_126("Pm-126", "Pm", 4, 65, 61, 126, 8008.0, "B-", std::numeric_limits<double>::quiet_NaN(), 125.957327, 537.0);
    static const Species Rh_127("Rh-127", "Rh", 37, 82, 45, 127, 8060.0, "B-", 13490.0, 126.963789, 644.0);
    static const Species Pd_127("Pd-127", "Pd", 35, 81, 46, 127, 8160.0, "B-", 11429.0, 126.949307, 537.0);
    static const Species Ag_127("Ag-127", "Ag", 33, 80, 47, 127, 8244.0, "B-", 10092.0, 126.937037, 215.0);
    static const Species Cd_127("Cd-127", "Cd", 31, 79, 48, 127, 8316.8971, "B-", 8138.6978, 126.926203291, 6.656);
    static const Species In_127("In-127", "In", 29, 78, 49, 127, 8374.8212, "B-", 6589.681, 126.91746604, 10.736);
    static const Species Sn_127("Sn-127", "Sn", 27, 77, 50, 127, 8420.5482, "B-", 3228.716, 126.910391726, 9.904);
    static const Species Sb_127("Sb-127", "Sb", 25, 76, 51, 127, 8439.811, "B-", 1582.203, 126.906925557, 5.457);
    static const Species Te_127("Te-127", "Te", 23, 75, 52, 127, 8446.109, "B-", 702.7199, 126.905226993, 1.465);
    static const Species I_127("I-127", "I", 21, 74, 53, 127, 8445.482, "B-", -662.3336, 126.904472592, 3.887);
    static const Species Xe_127("Xe-127", "Xe", 19, 73, 54, 127, 8434.1066, "B-", -2080.8562, 126.905183636, 4.388);
    static const Species Cs_127("Cs-127", "Cs", 17, 72, 55, 127, 8411.5617, "B-", -3422.0719, 126.907417527, 5.987);
    static const Species Ba_127("Ba-127", "Ba", 15, 71, 56, 127, 8378.456, "B-", -4921.8386, 126.911091272, 12.192);
    static const Species La_127("La-127", "La", 13, 70, 57, 127, 8333.5412, "B-", -5916.7727, 126.916375083, 27.912);
    static const Species Ce_127("Ce-127", "Ce", 11, 69, 58, 127, 8280.7922, "B-", -7436.0, 126.922727, 31.0);
    static const Species Pr_127("Pr-127", "Pr", 9, 68, 59, 127, 8216.0, "B-", -8633.0, 126.93071, 210.0);
    static const Species Nd_127("Nd-127", "Nd", 7, 67, 60, 127, 8142.0, "B-", -10600.0, 126.939978, 322.0);
    static const Species Pm_127("Pm-127", "Pm", 5, 66, 61, 127, 8052.0, "B-", std::numeric_limits<double>::quiet_NaN(), 126.951358, 429.0);
    static const Species Rh_128("Rh-128", "Rh", 38, 83, 45, 128, 8010.0, "B-", 17050.0, 127.970649, 322.0);
    static const Species Pd_128("Pd-128", "Pd", 36, 82, 46, 128, 8137.0, "B-", 10320.0, 127.952345, 537.0);
    static const Species Ag_128("Ag-128", "Ag", 34, 81, 47, 128, 8211.0, "B-", 12528.0, 127.941266, 322.0);
    static const Species Cd_128("Cd-128", "Cd", 32, 80, 48, 128, 8303.2367, "B-", 6951.8716, 127.927816778, 6.905);
    static const Species In_128("In-128", "In", 30, 79, 49, 128, 8351.4361, "B-", 9171.3131, 127.920353637, 1.419);
    static const Species Sn_128("Sn-128", "Sn", 28, 78, 50, 128, 8416.9749, "B-", 1268.4219, 127.910507828, 18.982);
    static const Species Sb_128("Sb-128", "Sb", 26, 77, 51, 128, 8420.7724, "B-", 4363.9419, 127.909146121, 20.169);
    static const Species Te_128("Te-128", "Te", 24, 76, 52, 128, 8448.7536, "B-", -1255.7634, 127.904461237, 0.758);
    static const Species I_128("I-128", "I", 22, 75, 53, 128, 8432.8309, "B-", 2122.5041, 127.905809355, 3.887);
    static const Species Xe_128("Xe-128", "Xe", 20, 74, 54, 128, 8443.3008, "B-", -3928.7617, 127.90353075341, 0.00558);
    static const Species Cs_128("Cs-128", "Cs", 18, 73, 55, 128, 8406.4953, "B-", -562.6171, 127.907748452, 5.771);
    static const Species Ba_128("Ba-128", "Ba", 16, 72, 56, 128, 8395.9878, "B-", -6743.7167, 127.908352446, 1.728);
    static const Species La_128("La-128", "La", 14, 71, 57, 128, 8337.1904, "B-", -3091.5136, 127.915592123, 58.452);
    static const Species Ce_128("Ce-128", "Ce", 12, 70, 58, 128, 8306.9259, "B-", -9203.1617, 127.918911, 30.0);
    static const Species Pr_128("Pr-128", "Pr", 10, 69, 59, 128, 8228.9141, "B-", -5800.0, 127.928791, 32.0);
    static const Species Nd_128("Nd-128", "Nd", 8, 68, 60, 128, 8177.0, "B-", -12311.0, 127.935018, 215.0);
    static const Species Pm_128("Pm-128", "Pm", 6, 67, 61, 128, 8075.0, "B-", -9070.0, 127.948234, 322.0);
    static const Species Sm_128("Sm-128", "Sm", 4, 66, 62, 128, 7998.0, "B-", std::numeric_limits<double>::quiet_NaN(), 127.957971, 537.0);
    static const Species Pd_129("Pd-129", "Pd", 37, 83, 46, 129, 8086.0, "B-", 13990.0, 128.959334, 644.0);
    static const Species Ag_129("Ag-129", "Ag", 35, 82, 47, 129, 8188.0, "B-", 11252.0, 128.944315, 429.0);
    static const Species Cd_129("Cd-129", "Cd", 33, 81, 48, 129, 8269.5311, "B-", 9712.7471, 128.932235597, 5.7);
    static const Species In_129("In-129", "In", 31, 80, 49, 129, 8338.759, "B-", 7755.7081, 128.921808534, 2.116);
    static const Species Sn_129("Sn-129", "Sn", 29, 79, 50, 129, 8392.8161, "B-", 4038.7874, 128.91348244, 18.54);
    static const Species Sb_129("Sb-129", "Sb", 27, 78, 51, 129, 8418.0598, "B-", 2375.5, 128.909146623, 22.786);
    static const Species Te_129("Te-129", "Te", 25, 77, 52, 129, 8430.4098, "B-", 1502.2919, 128.906596419, 0.763);
    static const Species I_129("I-129", "I", 23, 76, 53, 129, 8435.9908, "B-", 188.8936, 128.904983643, 3.385);
    static const Species Xe_129("Xe-129", "Xe", 21, 75, 54, 129, 8431.3904, "B-", -1197.0197, 128.90478085742, 0.00542);
    static const Species Cs_129("Cs-129", "Cs", 19, 74, 55, 129, 8416.0465, "B-", -2438.1843, 128.90606591, 4.888);
    static const Species Ba_129("Ba-129", "Ba", 17, 73, 56, 129, 8391.0811, "B-", -3737.3247, 128.908683409, 11.276);
    static const Species La_129("La-129", "La", 15, 72, 57, 129, 8356.0449, "B-", -5036.037, 128.912695592, 22.913);
    static const Species Ce_129("Ce-129", "Ce", 13, 71, 58, 129, 8310.9411, "B-", -6513.9383, 128.918102, 30.0);
    static const Species Pr_129("Pr-129", "Pr", 11, 70, 59, 129, 8254.3808, "B-", -7399.0, 128.925095, 32.0);
    static const Species Nd_129("Nd-129", "Nd", 9, 69, 60, 129, 8191.0, "B-", -9195.0, 128.933038, 217.0);
    static const Species Pm_129("Pm-129", "Pm", 7, 68, 61, 129, 8114.0, "B-", -10850.0, 128.942909, 322.0);
    static const Species Sm_129("Sm-129", "Sm", 5, 67, 62, 129, 8023.0, "B-", std::numeric_limits<double>::quiet_NaN(), 128.954557, 537.0);
    static const Species Pd_130("Pd-130", "Pd", 38, 84, 46, 130, 8046.0, "B-", 13168.0, 129.964863, 322.0);
    static const Species Ag_130("Ag-130", "Ag", 36, 83, 47, 130, 8142.0, "B-", 15220.0, 129.950727, 455.0);
    static const Species Cd_130("Cd-130", "Cd", 34, 82, 48, 130, 8252.5868, "B-", 8788.9322, 129.934387563, 24.0);
    static const Species In_130("In-130", "In", 32, 81, 49, 130, 8314.176, "B-", 10225.687, 129.924952257, 1.921);
    static const Species Sn_130("Sn-130", "Sn", 30, 80, 50, 130, 8386.817, "B-", 2153.4702, 129.913974531, 2.01);
    static const Species Sb_130("Sb-130", "Sb", 28, 79, 51, 130, 8397.3641, "B-", 5067.2728, 129.911662686, 15.257);
    static const Species Te_130("Te-130", "Te", 26, 78, 52, 130, 8430.3251, "B-", -416.7716, 129.906222745, 0.011);
    static const Species I_130("I-130", "I", 24, 77, 53, 130, 8421.1011, "B-", 2944.2864, 129.906670168, 3.385);
    static const Species Xe_130("Xe-130", "Xe", 22, 76, 54, 130, 8437.7314, "B-", -2980.7199, 129.903509346, 0.01);
    static const Species Cs_130("Cs-130", "Cs", 20, 75, 55, 130, 8408.7848, "B-", 357.0219, 129.906709281, 8.971);
    static const Species Ba_130("Ba-130", "Ba", 18, 74, 56, 130, 8405.513, "B-", -5629.4021, 129.906326002, 0.308);
    static const Species La_130("La-130", "La", 16, 73, 57, 130, 8356.1919, "B-", -2204.4611, 129.912369413, 27.854);
    static const Species Ce_130("Ce-130", "Ce", 14, 72, 58, 130, 8333.2164, "B-", -8247.4488, 129.914736, 30.0);
    static const Species Pr_130("Pr-130", "Pr", 12, 71, 59, 130, 8263.7565, "B-", -4579.225, 129.92359, 69.0);
    static const Species Nd_130("Nd-130", "Nd", 10, 70, 60, 130, 8222.5136, "B-", -11127.0, 129.928506, 30.0);
    static const Species Pm_130("Pm-130", "Pm", 8, 69, 61, 130, 8131.0, "B-", -7770.0, 129.940451, 215.0);
    static const Species Sm_130("Sm-130", "Sm", 6, 68, 62, 130, 8065.0, "B-", -14187.0, 129.948792, 429.0);
    static const Species Eu_130("Eu-130", "Eu", 4, 67, 63, 130, 7950.0, "B-", std::numeric_limits<double>::quiet_NaN(), 129.964022, 578.0);
    static const Species Pd_131("Pd-131", "Pd", 39, 85, 46, 131, 7993.0, "B-", 15010.0, 130.972367, 322.0);
    static const Species Ag_131("Ag-131", "Ag", 37, 84, 47, 131, 8102.0, "B-", 14462.0, 130.956253, 537.0);
    static const Species Cd_131("Cd-131", "Cd", 35, 83, 48, 131, 8206.1204, "B-", 12812.6089, 130.94072774, 20.653);
    static const Species In_131("In-131", "In", 33, 82, 49, 131, 8297.9544, "B-", 9240.2095, 130.926972839, 2.367);
    static const Species Sn_131("Sn-131", "Sn", 31, 81, 50, 131, 8362.5183, "B-", 4716.8328, 130.917053067, 3.887);
    static const Species Sb_131("Sb-131", "Sb", 29, 80, 51, 131, 8392.5525, "B-", 3229.6099, 130.911989339, 2.236);
    static const Species Te_131("Te-131", "Te", 27, 79, 52, 131, 8411.2339, "B-", 2231.7057, 130.90852221, 0.065);
    static const Species I_131("I-131", "I", 25, 78, 53, 131, 8422.2977, "B-", 970.8477, 130.906126375, 0.649);
    static const Species Xe_131("Xe-131", "Xe", 23, 77, 54, 131, 8423.7367, "B-", -358.0009, 130.90508412808, 0.00549);
    static const Species Cs_131("Cs-131", "Cs", 21, 76, 55, 131, 8415.0317, "B-", -1376.6158, 130.905468457, 0.19);
    static const Species Ba_131("Ba-131", "Ba", 19, 75, 56, 131, 8398.5511, "B-", -2909.6936, 130.906946315, 0.445);
    static const Species La_131("La-131", "La", 17, 74, 57, 131, 8370.3676, "B-", -4060.8167, 130.91007, 30.0);
    static const Species Ce_131("Ce-131", "Ce", 15, 73, 58, 131, 8333.3969, "B-", -5407.7842, 130.914429465, 35.214);
    static const Species Pr_131("Pr-131", "Pr", 13, 72, 59, 131, 8286.1439, "B-", -6532.6235, 130.92023496, 50.451);
    static const Species Nd_131("Nd-131", "Nd", 11, 71, 60, 131, 8230.3045, "B-", -7998.0, 130.92724802, 29.541);
    static const Species Pm_131("Pm-131", "Pm", 9, 70, 61, 131, 8163.0, "B-", -9490.0, 130.935834, 215.0);
    static const Species Sm_131("Sm-131", "Sm", 7, 69, 62, 131, 8085.0, "B-", -10816.0, 130.946022, 429.0);
    static const Species Eu_131("Eu-131", "Eu", 5, 68, 63, 131, 7996.0, "B-", std::numeric_limits<double>::quiet_NaN(), 130.957634, 429.0);
    static const Species Ag_132("Ag-132", "Ag", 38, 85, 47, 132, 8053.0, "B-", 16065.0, 131.96307, 537.0);
    static const Species Cd_132("Cd-132", "Cd", 36, 84, 48, 132, 8169.1421, "B-", 11946.1243, 131.945823136, 64.485);
    static const Species In_132("In-132", "In", 34, 83, 49, 132, 8253.7162, "B-", 14135.0, 131.932998444, 64.447);
    static const Species Sn_132("Sn-132", "Sn", 32, 82, 50, 132, 8354.8726, "B-", 3088.728, 131.917823898, 2.121);
    static const Species Sb_132("Sb-132", "Sb", 30, 81, 51, 132, 8372.3452, "B-", 5552.9155, 131.914508013, 2.648);
    static const Species Te_132("Te-132", "Te", 28, 80, 52, 132, 8408.4859, "B-", 515.3046, 131.908546713, 3.742);
    static const Species I_132("I-132", "I", 26, 79, 53, 132, 8406.4628, "B-", 3575.4729, 131.907993511, 4.364);
    static const Species Xe_132("Xe-132", "Xe", 24, 78, 54, 132, 8427.6229, "B-", -2126.2813, 131.90415508346, 0.00544);
    static const Species Cs_132("Cs-132", "Cs", 22, 77, 55, 132, 8405.5878, "B-", 1282.2099, 131.90643774, 1.112);
    static const Species Ba_132("Ba-132", "Ba", 20, 76, 56, 132, 8409.3747, "B-", -4711.3256, 131.905061231, 1.13);
    static const Species La_132("La-132", "La", 18, 75, 57, 132, 8367.7559, "B-", -1254.8898, 131.910119047, 39.032);
    static const Species Ce_132("Ce-132", "Ce", 16, 74, 58, 132, 8352.3223, "B-", -7241.224, 131.911466226, 21.907);
    static const Species Pr_132("Pr-132", "Pr", 14, 73, 59, 132, 8291.5377, "B-", -3801.6487, 131.91924, 31.0);
    static const Species Nd_132("Nd-132", "Nd", 12, 72, 60, 132, 8256.8104, "B-", -9798.0, 131.923321237, 25.985);
    static const Species Pm_132("Pm-132", "Pm", 10, 71, 61, 132, 8177.0, "B-", -6488.0, 131.93384, 160.0);
    static const Species Sm_132("Sm-132", "Sm", 8, 70, 62, 132, 8122.0, "B-", -12939.0, 131.940805, 322.0);
    static const Species Eu_132("Eu-132", "Eu", 6, 69, 63, 132, 8018.0, "B-", std::numeric_limits<double>::quiet_NaN(), 131.954696, 429.0);
    static const Species Ag_133("Ag-133", "Ag", 39, 86, 47, 133, 8013.0, "B-", 15059.0, 132.968781, 537.0);
    static const Species Cd_133("Cd-133", "Cd", 37, 85, 48, 133, 8121.0, "B-", 13550.0, 132.952614, 215.0);
    static const Species In_133("In-133", "In", 35, 84, 49, 133, 8217.0, "B-", 13184.0, 132.938067, 215.0);
    static const Species Sn_133("Sn-133", "Sn", 33, 83, 50, 133, 8310.089, "B-", 8049.6228, 132.923913753, 2.043);
    static const Species Sb_133("Sb-133", "Sb", 31, 82, 51, 133, 8364.7302, "B-", 4013.6198, 132.915272128, 3.357);
    static const Species Te_133("Te-133", "Te", 29, 81, 52, 133, 8389.0255, "B-", 2920.169, 132.91096333, 2.218);
    static const Species I_133("I-133", "I", 27, 80, 53, 133, 8405.0993, "B-", 1786.2812, 132.9078284, 6.335);
    static const Species Xe_133("Xe-133", "Xe", 25, 79, 54, 133, 8412.6477, "B-", 427.36, 132.905910748, 2.576);
    static const Species Cs_133("Cs-133", "Cs", 23, 78, 55, 133, 8409.9786, "B-", -517.431, 132.905451958, 0.008);
    static const Species Ba_133("Ba-133", "Ba", 21, 77, 56, 133, 8400.2059, "B-", -2059.1203, 132.906007443, 1.065);
    static const Species La_133("La-133", "La", 19, 76, 57, 133, 8378.8415, "B-", -3076.1685, 132.908218, 30.0);
    static const Species Ce_133("Ce-133", "Ce", 17, 75, 58, 133, 8349.8301, "B-", -4480.6319, 132.911520402, 17.557);
    static const Species Pr_133("Pr-133", "Pr", 15, 74, 59, 133, 8310.2588, "B-", -5605.2112, 132.916330558, 13.416);
    static const Species Nd_133("Nd-133", "Nd", 13, 73, 60, 133, 8262.232, "B-", -6924.7272, 132.922348, 50.0);
    static const Species Pm_133("Pm-133", "Pm", 11, 72, 61, 133, 8204.2841, "B-", -8177.0, 132.929782, 54.0);
    static const Species Sm_133("Sm-133", "Sm", 9, 71, 62, 133, 8137.0, "B-", -9995.0, 132.93856, 320.0);
    static const Species Eu_133("Eu-133", "Eu", 7, 70, 63, 133, 8056.0, "B-", -11176.0, 132.94929, 320.0);
    static const Species Gd_133("Gd-133", "Gd", 5, 69, 64, 133, 7966.0, "B-", std::numeric_limits<double>::quiet_NaN(), 132.961288, 537.0);
    static const Species Cd_134("Cd-134", "Cd", 38, 86, 48, 134, 8086.0, "B-", 12510.0, 133.957638, 322.0);
    static const Species In_134("In-134", "In", 36, 85, 49, 134, 8173.0, "B-", 14464.0, 133.944208, 215.0);
    static const Species Sn_134("Sn-134", "Sn", 34, 84, 50, 134, 8275.1719, "B-", 7585.2453, 133.92868043, 3.4);
    static const Species Sb_134("Sb-134", "Sb", 32, 83, 51, 134, 8325.9398, "B-", 8514.7483, 133.920537334, 3.3);
    static const Species Te_134("Te-134", "Te", 30, 82, 52, 134, 8383.6442, "B-", 1509.6875, 133.911396376, 2.948);
    static const Species I_134("I-134", "I", 28, 81, 53, 134, 8389.0722, "B-", 4082.3946, 133.90977566, 5.213);
    static const Species Xe_134("Xe-134", "Xe", 26, 80, 54, 134, 8413.6994, "B-", -1234.6691, 133.90539303, 0.006);
    static const Species Cs_134("Cs-134", "Cs", 24, 79, 55, 134, 8398.647, "B-", 2058.8368, 133.906718501, 0.017);
    static const Species Ba_134("Ba-134", "Ba", 22, 78, 56, 134, 8408.1731, "B-", -3731.3434, 133.904508249, 0.269);
    static const Species La_134("La-134", "La", 20, 77, 57, 134, 8374.4888, "B-", -385.7605, 133.908514011, 21.395);
    static const Species Ce_134("Ce-134", "Ce", 18, 76, 58, 134, 8365.7716, "B-", -6304.8987, 133.908928142, 21.886);
    static const Species Pr_134("Pr-134", "Pr", 16, 75, 59, 134, 8312.8817, "B-", -2881.5569, 133.915696729, 21.81);
    static const Species Nd_134("Nd-134", "Nd", 14, 74, 60, 134, 8285.5391, "B-", -8882.5343, 133.918790207, 12.686);
    static const Species Pm_134("Pm-134", "Pm", 12, 73, 61, 134, 8213.4131, "B-", -5388.0, 133.928326, 45.0);
    static const Species Sm_134("Sm-134", "Sm", 10, 72, 62, 134, 8167.0, "B-", -11576.0, 133.93411, 210.0);
    static const Species Eu_134("Eu-134", "Eu", 8, 71, 63, 134, 8075.0, "B-", -8271.0, 133.946537, 322.0);
    static const Species Gd_134("Gd-134", "Gd", 6, 70, 64, 134, 8008.0, "B-", std::numeric_limits<double>::quiet_NaN(), 133.955416, 429.0);
    static const Species Cd_135("Cd-135", "Cd", 39, 87, 48, 135, 8036.0, "B-", 14290.0, 134.964766, 429.0);
    static const Species In_135("In-135", "In", 37, 86, 49, 135, 8136.0, "B-", 13522.0, 134.949425, 322.0);
    static const Species Sn_135("Sn-135", "Sn", 35, 85, 50, 135, 8230.6877, "B-", 9058.08, 134.934908603, 3.3);
    static const Species Sb_135("Sb-135", "Sb", 33, 84, 51, 135, 8291.9894, "B-", 8038.4581, 134.925184354, 2.834);
    static const Species Te_135("Te-135", "Te", 31, 83, 52, 135, 8345.7384, "B-", 6050.3894, 134.916554715, 1.848);
    static const Species I_135("I-135", "I", 29, 82, 53, 135, 8384.7609, "B-", 2634.1851, 134.910059355, 2.211);
    static const Species Xe_135("Xe-135", "Xe", 27, 81, 54, 135, 8398.4783, "B-", 1168.5917, 134.907231441, 3.938);
    static const Species Cs_135("Cs-135", "Cs", 25, 80, 55, 135, 8401.3393, "B-", 268.6983, 134.905976907, 0.39);
    static const Species Ba_135("Ba-135", "Ba", 23, 79, 56, 135, 8397.5345, "B-", -1207.1973, 134.905688447, 0.263);
    static const Species La_135("La-135", "La", 21, 78, 57, 135, 8382.7972, "B-", -2027.1499, 134.906984427, 10.126);
    static const Species Ce_135("Ce-135", "Ce", 19, 77, 58, 135, 8361.9861, "B-", -3680.4357, 134.909160662, 11.021);
    static const Species Pr_135("Pr-135", "Pr", 17, 76, 59, 135, 8328.9284, "B-", -4722.2522, 134.913111772, 12.686);
    static const Species Nd_135("Nd-135", "Nd", 15, 75, 60, 135, 8288.1536, "B-", -6151.2907, 134.918181318, 20.534);
    static const Species Pm_135("Pm-135", "Pm", 13, 74, 61, 135, 8236.7933, "B-", -7205.1069, 134.92478499999999, 89.0);
    static const Species Sm_135("Sm-135", "Sm", 11, 73, 62, 135, 8177.627, "B-", -8709.0, 134.93252, 166.0);
    static const Species Eu_135("Eu-135", "Eu", 9, 72, 63, 135, 8107.0, "B-", -9898.0, 134.94187, 210.0);
    static const Species Gd_135("Gd-135", "Gd", 7, 71, 64, 135, 8028.0, "B-", -11197.0, 134.952496, 429.0);
    static const Species Tb_135("Tb-135", "Tb", 5, 70, 65, 135, 7939.0, "B-", std::numeric_limits<double>::quiet_NaN(), 134.964516, 429.0);
    static const Species In_136("In-136", "In", 38, 87, 49, 136, 8091.0, "B-", 15200.0, 135.956017, 322.0);
    static const Species Sn_136("Sn-136", "Sn", 36, 86, 50, 136, 8197.0, "B-", 8337.0, 135.939699, 215.0);
    static const Species Sb_136("Sb-136", "Sb", 34, 85, 51, 136, 8252.2533, "B-", 9918.3897, 135.930749009, 6.258);
    static const Species Te_136("Te-136", "Te", 32, 84, 52, 136, 8319.4301, "B-", 5119.9455, 135.92010118, 2.448);
    static const Species I_136("I-136", "I", 30, 83, 53, 136, 8351.3242, "B-", 6883.9455, 135.914604693, 15.231);
    static const Species Xe_136("Xe-136", "Xe", 28, 82, 54, 136, 8396.1889, "B-", -90.3151, 135.907214474, 0.007);
    static const Species Cs_136("Cs-136", "Cs", 26, 81, 55, 136, 8389.7723, "B-", 2548.2241, 135.907311431, 2.01);
    static const Species Ba_136("Ba-136", "Ba", 24, 80, 56, 136, 8402.7566, "B-", -2849.5915, 135.9045758, 0.262);
    static const Species La_136("La-136", "La", 22, 79, 57, 136, 8376.0512, "B-", 471.0621, 135.907634962, 57.081);
    static const Species Ce_136("Ce-136", "Ce", 20, 78, 58, 136, 8373.7624, "B-", -5168.129, 135.907129256, 0.348);
    static const Species Pr_136("Pr-136", "Pr", 18, 77, 59, 136, 8330.0089, "B-", -2141.1234, 135.91267747, 12.296);
    static const Species Nd_136("Nd-136", "Nd", 16, 76, 60, 136, 8308.5127, "B-", -8029.3747, 135.914976061, 12.686);
    static const Species Pm_136("Pm-136", "Pm", 14, 75, 61, 136, 8243.7207, "B-", -4359.0237, 135.923595949, 74.152);
    static const Species Sm_136("Sm-136", "Sm", 12, 74, 62, 136, 8205.9165, "B-", -10567.0, 135.928275553, 13.416);
    static const Species Eu_136("Eu-136", "Eu", 10, 73, 63, 136, 8122.0, "B-", -7154.0, 135.93962, 210.0);
    static const Species Gd_136("Gd-136", "Gd", 8, 72, 64, 136, 8064.0, "B-", -13190.0, 135.9473, 320.0);
    static const Species Tb_136("Tb-136", "Tb", 6, 71, 65, 136, 7961.0, "B-", std::numeric_limits<double>::quiet_NaN(), 135.96146, 537.0);
    static const Species In_137("In-137", "In", 39, 88, 49, 137, 8053.0, "B-", 14320.0, 136.961535, 429.0);
    static const Species Sn_137("Sn-137", "Sn", 37, 87, 50, 137, 8152.0, "B-", 9911.0, 136.946162, 322.0);
    static const Species Sb_137("Sb-137", "Sb", 35, 86, 51, 137, 8218.4764, "B-", 9243.3698, 136.935522519, 56.0);
    static const Species Te_137("Te-137", "Te", 33, 85, 52, 137, 8280.2357, "B-", 7052.5063, 136.925599354, 2.254);
    static const Species I_137("I-137", "I", 31, 84, 53, 137, 8326.0033, "B-", 6027.1455, 136.918028178, 9.0);
    static const Species Xe_137("Xe-137", "Xe", 29, 83, 54, 137, 8364.2865, "B-", 4162.3582, 136.911557771, 0.111);
    static const Species Cs_137("Cs-137", "Cs", 27, 82, 55, 137, 8388.9581, "B-", 1175.6285, 136.907089296, 0.324);
    static const Species Ba_137("Ba-137", "Ba", 25, 81, 56, 137, 8391.8288, "B-", -580.5356, 136.905827207, 0.266);
    static const Species La_137("La-137", "La", 23, 80, 57, 137, 8381.8807, "B-", -1222.1, 136.906450438, 1.76);
    static const Species Ce_137("Ce-137", "Ce", 21, 79, 58, 137, 8367.2497, "B-", -2716.9515, 136.907762416, 0.386);
    static const Species Pr_137("Pr-137", "Pr", 19, 78, 59, 137, 8341.7074, "B-", -3617.8447, 136.910679183, 8.733);
    static const Species Nd_137("Nd-137", "Nd", 17, 77, 60, 137, 8309.5892, "B-", -5511.1108, 136.914563099, 12.586);
    static const Species Pm_137("Pm-137", "Pm", 15, 76, 61, 137, 8263.6516, "B-", -6081.2029, 136.920479519, 14.0);
    static const Species Sm_137("Sm-137", "Sm", 13, 75, 62, 137, 8213.5527, "B-", -7845.7518, 136.927007959, 30.718);
    static const Species Eu_137("Eu-137", "Eu", 11, 74, 63, 137, 8150.5738, "B-", -8932.0, 136.935430719, 4.7);
    static const Species Gd_137("Gd-137", "Gd", 9, 73, 64, 137, 8080.0, "B-", -10246.0, 136.94502, 320.0);
    static const Species Tb_137("Tb-137", "Tb", 7, 72, 65, 137, 7999.0, "B-", std::numeric_limits<double>::quiet_NaN(), 136.95602, 430.0);
    static const Species Sn_138("Sn-138", "Sn", 38, 88, 50, 138, 8118.0, "B-", 9140.0, 137.951143, 429.0);
    static const Species Sb_138("Sb-138", "Sb", 36, 87, 51, 138, 8178.0, "B-", 11046.0, 137.941331, 322.0);
    static const Species Te_138("Te-138", "Te", 34, 86, 52, 138, 8252.5786, "B-", 6283.9149, 137.929472452, 4.065);
    static const Species I_138("I-138", "I", 32, 85, 53, 138, 8292.445, "B-", 7992.3346, 137.922726392, 6.4);
    static const Species Xe_138("Xe-138", "Xe", 30, 84, 54, 138, 8344.6913, "B-", 2914.7839, 137.914146268, 3.01);
    static const Species Cs_138("Cs-138", "Cs", 28, 83, 55, 138, 8360.1437, "B-", 5374.7776, 137.911017119, 9.831);
    static const Species Ba_138("Ba-138", "Ba", 26, 82, 56, 138, 8393.4222, "B-", -1748.3977, 137.905247059, 0.267);
    static const Species La_138("La-138", "La", 24, 81, 57, 138, 8375.0835, "B-", 1052.4585, 137.907124041, 0.446);
    static const Species Ce_138("Ce-138", "Ce", 22, 80, 58, 138, 8377.0408, "B-", -4437.0, 137.90599418, 0.536);
    static const Species Pr_138("Pr-138", "Pr", 20, 79, 59, 138, 8339.2195, "B-", -1111.6847, 137.910757495, 10.748);
    static const Species Nd_138("Nd-138", "Nd", 18, 78, 60, 138, 8325.4946, "B-", -7102.8119, 137.911950938, 12.456);
    static const Species Pm_138("Pm-138", "Pm", 16, 77, 61, 138, 8268.3558, "B-", -3416.5976, 137.919576119, 12.456);
    static const Species Sm_138("Sm-138", "Sm", 14, 76, 62, 138, 8237.9286, "B-", -9748.0968, 137.923243988, 12.686);
    static const Species Eu_138("Eu-138", "Eu", 12, 75, 63, 138, 8161.6211, "B-", -6090.0, 137.933709, 30.0);
    static const Species Gd_138("Gd-138", "Gd", 10, 74, 64, 138, 8112.0, "B-", -12059.0, 137.940247, 215.0);
    static const Species Tb_138("Tb-138", "Tb", 8, 73, 65, 138, 8019.0, "B-", -8669.0, 137.953193, 322.0);
    static const Species Dy_138("Dy-138", "Dy", 6, 72, 66, 138, 7950.0, "B-", std::numeric_limits<double>::quiet_NaN(), 137.9625, 540.0);
    static const Species Sn_139("Sn-139", "Sn", 39, 89, 50, 139, 8073.0, "B-", 10740.0, 138.957799, 429.0);
    static const Species Sb_139("Sb-139", "Sb", 37, 88, 51, 139, 8144.0, "B-", 10155.0, 138.946269, 429.0);
    static const Species Te_139("Te-139", "Te", 35, 87, 52, 139, 8211.7716, "B-", 8265.8835, 138.935367191, 3.8);
    static const Species I_139("I-139", "I", 33, 86, 53, 139, 8265.61, "B-", 7173.6224, 138.9264934, 4.3);
    static const Species Xe_139("Xe-139", "Xe", 31, 85, 54, 139, 8311.5904, "B-", 5056.5023, 138.9187922, 2.3);
    static const Species Cs_139("Cs-139", "Cs", 29, 84, 55, 139, 8342.3397, "B-", 4212.8293, 138.913363822, 3.364);
    static const Species Ba_139("Ba-139", "Ba", 27, 83, 56, 139, 8367.0194, "B-", 2308.4632, 138.908841164, 0.271);
    static const Species La_139("La-139", "La", 25, 82, 57, 139, 8377.9987, "B-", -264.6396, 138.906362927, 0.651);
    static const Species Ce_139("Ce-139", "Ce", 23, 81, 58, 139, 8370.4664, "B-", -2129.089, 138.906647029, 2.242);
    static const Species Pr_139("Pr-139", "Pr", 21, 80, 59, 139, 8349.5208, "B-", -2811.7226, 138.9089327, 3.917);
    static const Species Nd_139("Nd-139", "Nd", 19, 79, 60, 139, 8323.6642, "B-", -4515.9014, 138.911951208, 29.545);
    static const Species Pm_139("Pm-139", "Pm", 17, 78, 61, 139, 8285.5473, "B-", -5120.7993, 138.916799228, 14.587);
    static const Species Sm_139("Sm-139", "Sm", 15, 77, 62, 139, 8243.0786, "B-", -6982.1777, 138.922296631, 11.684);
    static const Species Eu_139("Eu-139", "Eu", 13, 76, 63, 139, 8187.2187, "B-", -7767.0, 138.929792307, 14.117);
    static const Species Gd_139("Gd-139", "Gd", 11, 75, 64, 139, 8126.0, "B-", -9501.0, 138.93813, 210.0);
    static const Species Tb_139("Tb-139", "Tb", 9, 74, 65, 139, 8052.0, "B-", -10430.0, 138.94833, 320.0);
    static const Species Dy_139("Dy-139", "Dy", 7, 73, 66, 139, 7971.0, "B-", std::numeric_limits<double>::quiet_NaN(), 138.959527, 537.0);
    static const Species Sn_140("Sn-140", "Sn", 40, 90, 50, 140, 8038.0, "B-", 9900.0, 139.962973, 322.0);
    static const Species Sb_140("Sb-140", "Sb", 38, 89, 51, 140, 8103.0, "B-", 11977.0, 139.952345, 644.0);
    static const Species Te_140("Te-140", "Te", 36, 88, 52, 140, 8183.3567, "B-", 7238.7734, 139.939487057, 15.434);
    static const Species I_140("I-140", "I", 34, 87, 53, 140, 8229.474, "B-", 9380.2388, 139.931715914, 13.0);
    static const Species Xe_140("Xe-140", "Xe", 32, 86, 54, 140, 8290.8875, "B-", 4063.2768, 139.921645814, 2.5);
    static const Species Cs_140("Cs-140", "Cs", 30, 85, 55, 140, 8314.3227, "B-", 6218.1669, 139.917283707, 8.801);
    static const Species Ba_140("Ba-140", "Ba", 28, 84, 56, 140, 8353.15, "B-", 1044.1542, 139.910608231, 8.48);
    static const Species La_140("La-140", "La", 26, 83, 57, 140, 8355.0201, "B-", 3762.1676, 139.909487285, 0.651);
    static const Species Ce_140("Ce-140", "Ce", 24, 82, 58, 140, 8376.3045, "B-", -3388.0, 139.905448433, 1.409);
    static const Species Pr_140("Pr-140", "Pr", 22, 81, 59, 140, 8346.5163, "B-", -428.9806, 139.9090856, 6.593);
    static const Species Nd_140("Nd-140", "Nd", 20, 80, 60, 140, 8337.864, "B-", -6045.2, 139.90954613, 3.5);
    static const Species Pm_140("Pm-140", "Pm", 18, 79, 61, 140, 8289.0958, "B-", -2756.101, 139.916035918, 26.001);
    static const Species Sm_140("Sm-140", "Sm", 16, 78, 62, 140, 8263.8211, "B-", -8470.0, 139.918994714, 13.416);
    static const Species Eu_140("Eu-140", "Eu", 14, 77, 63, 140, 8197.733, "B-", -5203.6675, 139.928087633, 55.328);
    static const Species Gd_140("Gd-140", "Gd", 12, 76, 64, 140, 8154.9757, "B-", -11300.0, 139.933674, 30.0);
    static const Species Tb_140("Tb-140", "Tb", 10, 75, 65, 140, 8068.6732, "B-", -7652.0, 139.945805048, 859.359);
    static const Species Dy_140("Dy-140", "Dy", 8, 74, 66, 140, 8008.0, "B-", -13513.0, 139.95402, 430.0);
    static const Species Ho_140("Ho-140", "Ho", 6, 73, 67, 140, 7906.0, "B-", std::numeric_limits<double>::quiet_NaN(), 139.968526, 537.0);
    static const Species Sb_141("Sb-141", "Sb", 39, 90, 51, 141, 8069.0, "B-", 11129.0, 140.957552, 537.0);
    static const Species Te_141("Te-141", "Te", 37, 89, 52, 141, 8142.0, "B-", 9257.0, 140.945604, 429.0);
    static const Species I_141("I-141", "I", 35, 88, 53, 141, 8202.2562, "B-", 8270.643, 140.935666081, 17.0);
    static const Species Xe_141("Xe-141", "Xe", 33, 87, 54, 141, 8255.3647, "B-", 6280.0423, 140.926787181, 3.1);
    static const Species Cs_141("Cs-141", "Cs", 31, 86, 55, 141, 8294.3554, "B-", 5255.141, 140.920045279, 9.871);
    static const Species Ba_141("Ba-141", "Ba", 29, 85, 56, 141, 8326.0774, "B-", 3197.3515, 140.914403653, 5.709);
    static const Species La_141("La-141", "La", 27, 84, 57, 141, 8343.205, "B-", 2501.2141, 140.910971155, 4.43);
    static const Species Ce_141("Ce-141", "Ce", 25, 83, 58, 141, 8355.3956, "B-", 583.4758, 140.908285991, 1.411);
    static const Species Pr_141("Pr-141", "Pr", 23, 82, 59, 141, 8353.9852, "B-", -1823.0137, 140.907659604, 1.607);
    static const Species Nd_141("Nd-141", "Nd", 21, 81, 60, 141, 8335.5074, "B-", -3668.5879, 140.90961669, 3.417);
    static const Species Pm_141("Pm-141", "Pm", 19, 80, 61, 141, 8303.9405, "B-", -4588.9724, 140.913555081, 15.0);
    static const Species Sm_141("Sm-141", "Sm", 17, 79, 62, 141, 8265.846, "B-", -6008.3127, 140.918481545, 9.162);
    static const Species Eu_141("Eu-141", "Eu", 15, 78, 63, 141, 8217.6853, "B-", -6701.4161, 140.924931734, 13.568);
    static const Species Gd_141("Gd-141", "Gd", 13, 77, 64, 141, 8164.609, "B-", -8683.388, 140.932126, 21.213);
    static const Species Tb_141("Tb-141", "Tb", 11, 76, 65, 141, 8097.4761, "B-", -9158.0, 140.941448, 113.0);
    static const Species Dy_141("Dy-141", "Dy", 9, 75, 66, 141, 8027.0, "B-", -11018.0, 140.95128, 320.0);
    static const Species Ho_141("Ho-141", "Ho", 7, 74, 67, 141, 7943.0, "B-", std::numeric_limits<double>::quiet_NaN(), 140.963108, 430.0);
    static const Species Sb_142("Sb-142", "Sb", 40, 91, 51, 142, 8027.0, "B-", 12939.0, 141.963918, 322.0);
    static const Species Te_142("Te-142", "Te", 38, 90, 52, 142, 8113.0, "B-", 8253.0, 141.950027, 537.0);
    static const Species I_142("I-142", "I", 36, 89, 53, 142, 8165.2517, "B-", 10426.6792, 141.941166595, 5.3);
    static const Species Xe_142("Xe-142", "Xe", 34, 88, 54, 142, 8233.1695, "B-", 5284.9078, 141.929973095, 2.9);
    static const Species Cs_142("Cs-142", "Cs", 32, 87, 55, 142, 8264.8777, "B-", 7327.7007, 141.924299514, 7.586);
    static const Species Ba_142("Ba-142", "Ba", 30, 86, 56, 142, 8310.9718, "B-", 2181.6932, 141.916432904, 6.355);
    static const Species La_142("La-142", "La", 28, 85, 57, 142, 8320.8263, "B-", 4508.9455, 141.91409076, 6.748);
    static const Species Ce_142("Ce-142", "Ce", 26, 84, 58, 142, 8347.07, "B-", -746.5288, 141.909250208, 2.623);
    static const Species Pr_142("Pr-142", "Pr", 24, 83, 59, 142, 8336.3032, "B-", 2163.6885, 141.91005164, 1.607);
    static const Species Nd_142("Nd-142", "Nd", 22, 82, 60, 142, 8346.031, "B-", -4808.5196, 141.907728824, 1.348);
    static const Species Pm_142("Pm-142", "Pm", 20, 81, 61, 142, 8306.6587, "B-", -2159.6062, 141.912890982, 25.33);
    static const Species Sm_142("Sm-142", "Sm", 18, 80, 62, 142, 8285.9407, "B-", -7673.0, 141.915209415, 2.002);
    static const Species Eu_142("Eu-142", "Eu", 16, 79, 63, 142, 8226.396, "B-", -4349.4074, 141.923446719, 32.268);
    static const Species Gd_142("Gd-142", "Gd", 14, 78, 64, 142, 8190.2569, "B-", -10400.0, 141.928116, 30.0);
    static const Species Tb_142("Tb-142", "Tb", 12, 77, 65, 142, 8111.508, "B-", -6440.0, 141.939280858, 752.079);
    static const Species Dy_142("Dy-142", "Dy", 10, 76, 66, 142, 8061.0, "B-", -12869.0, 141.946194, 782.0);
    static const Species Ho_142("Ho-142", "Ho", 8, 75, 67, 142, 7965.0, "B-", -9321.0, 141.96001, 430.0);
    static const Species Er_142("Er-142", "Er", 6, 74, 68, 142, 7893.0, "B-", std::numeric_limits<double>::quiet_NaN(), 141.970016, 537.0);
    static const Species Te_143("Te-143", "Te", 39, 91, 52, 143, 8070.0, "B-", 10259.0, 142.956489, 537.0);
    static const Species I_143("I-143", "I", 37, 90, 53, 143, 8137.0, "B-", 9413.0, 142.945475, 215.0);
    static const Species Xe_143("Xe-143", "Xe", 35, 89, 54, 143, 8196.8855, "B-", 7472.6365, 142.93536955, 5.0);
    static const Species Cs_143("Cs-143", "Cs", 33, 88, 55, 143, 8243.6707, "B-", 6261.6865, 142.927347346, 8.13);
    static const Species Ba_143("Ba-143", "Ba", 31, 87, 56, 143, 8281.9878, "B-", 4234.2623, 142.920625149, 7.253);
    static const Species La_143("La-143", "La", 29, 86, 57, 143, 8306.1271, "B-", 3434.9108, 142.916079482, 7.868);
    static const Species Ce_143("Ce-143", "Ce", 27, 85, 58, 143, 8324.6765, "B-", 1461.8214, 142.912391953, 2.621);
    static const Species Pr_143("Pr-143", "Pr", 25, 84, 59, 143, 8329.428, "B-", 934.1107, 142.910822624, 1.949);
    static const Species Nd_143("Nd-143", "Nd", 23, 83, 60, 143, 8330.4893, "B-", -1041.6463, 142.909819815, 1.347);
    static const Species Pm_143("Pm-143", "Pm", 21, 82, 61, 143, 8317.7341, "B-", -3443.5291, 142.910938068, 3.16);
    static const Species Sm_143("Sm-143", "Sm", 19, 81, 62, 143, 8288.1825, "B-", -5275.824, 142.914634848, 2.951);
    static const Species Eu_143("Eu-143", "Eu", 17, 80, 63, 143, 8245.8177, "B-", -6010.0, 142.920298678, 11.793);
    static const Species Gd_143("Gd-143", "Gd", 15, 79, 64, 143, 8198.3188, "B-", -7812.1185, 142.926750678, 215.032);
    static const Species Tb_143("Tb-143", "Tb", 13, 78, 65, 143, 8138.2176, "B-", -8250.2433, 142.935137332, 55.0);
    static const Species Dy_143("Dy-143", "Dy", 11, 77, 66, 143, 8075.0527, "B-", -10121.0, 142.943994332, 14.0);
    static const Species Ho_143("Ho-143", "Ho", 9, 76, 67, 143, 7999.0, "B-", -10887.0, 142.95486, 320.0);
    static const Species Er_143("Er-143", "Er", 7, 75, 68, 143, 7917.0, "B-", std::numeric_limits<double>::quiet_NaN(), 142.966548, 429.0);
    static const Species Te_144("Te-144", "Te", 40, 92, 52, 144, 8040.0, "B-", 9110.0, 143.961116, 322.0);
    static const Species I_144("I-144", "I", 38, 91, 53, 144, 8098.0, "B-", 11542.0, 143.951336, 429.0);
    static const Species Xe_144("Xe-144", "Xe", 36, 90, 54, 144, 8172.8845, "B-", 6399.0606, 143.938945076, 5.7);
    static const Species Cs_144("Cs-144", "Cs", 34, 89, 55, 144, 8211.8894, "B-", 8495.7679, 143.932075402, 21.612);
    static const Species Ba_144("Ba-144", "Ba", 32, 88, 56, 144, 8265.4549, "B-", 3082.53, 143.922954821, 7.661);
    static const Species La_144("La-144", "La", 30, 87, 57, 144, 8281.4283, "B-", 5582.2823, 143.919645589, 13.888);
    static const Species Ce_144("Ce-144", "Ce", 28, 86, 58, 144, 8314.7612, "B-", 318.6462, 143.913652763, 3.041);
    static const Species Pr_144("Pr-144", "Pr", 26, 85, 59, 144, 8311.5411, "B-", 2997.44, 143.913310682, 2.907);
    static const Species Nd_144("Nd-144", "Nd", 24, 84, 60, 144, 8326.9237, "B-", -2331.9117, 143.910092798, 1.346);
    static const Species Pm_144("Pm-144", "Pm", 22, 83, 61, 144, 8305.2969, "B-", 549.5096, 143.912596208, 3.126);
    static const Species Sm_144("Sm-144", "Sm", 20, 82, 62, 144, 8303.68, "B-", -6346.4515, 143.912006285, 1.566);
    static const Species Eu_144("Eu-144", "Eu", 18, 81, 63, 144, 8254.1744, "B-", -3859.6633, 143.918819481, 11.58);
    static const Species Gd_144("Gd-144", "Gd", 16, 80, 64, 144, 8221.9382, "B-", -9391.3235, 143.922963, 30.0);
    static const Species Tb_144("Tb-144", "Tb", 14, 79, 65, 144, 8151.2877, "B-", -5798.0965, 143.933045, 30.0);
    static const Species Dy_144("Dy-144", "Dy", 12, 78, 66, 144, 8105.5902, "B-", -11960.5706, 143.939269512, 7.7);
    static const Species Ho_144("Ho-144", "Ho", 10, 77, 67, 144, 8017.0977, "B-", -8002.0, 143.952109712, 9.1);
    static const Species Er_144("Er-144", "Er", 8, 76, 68, 144, 7956.0, "B-", -14448.0, 143.9607, 210.0);
    static const Species Tm_144("Tm-144", "Tm", 6, 75, 69, 144, 7850.0, "B-", std::numeric_limits<double>::quiet_NaN(), 143.976211, 429.0);
    static const Species Te_145("Te-145", "Te", 41, 93, 52, 145, 7998.0, "B-", 11120.0, 144.967783, 322.0);
    static const Species I_145("I-145", "I", 39, 92, 53, 145, 8069.0, "B-", 10363.0, 144.955845, 537.0);
    static const Species Xe_145("Xe-145", "Xe", 37, 91, 54, 145, 8135.0877, "B-", 8561.0867, 144.944719631, 12.0);
    static const Species Cs_145("Cs-145", "Cs", 35, 90, 55, 145, 8188.7342, "B-", 7461.7591, 144.935528927, 9.733);
    static const Species Ba_145("Ba-145", "Ba", 33, 89, 56, 145, 8234.7991, "B-", 5319.1425, 144.9275184, 9.1);
    static const Species La_145("La-145", "La", 31, 88, 57, 145, 8266.0873, "B-", 4231.7331, 144.921808065, 13.17);
    static const Species Ce_145("Ce-145", "Ce", 29, 87, 58, 145, 8289.8762, "B-", 2558.9318, 144.917265113, 36.393);
    static const Species Pr_145("Pr-145", "Pr", 27, 86, 59, 145, 8302.1285, "B-", 1806.014, 144.914517987, 7.674);
    static const Species Nd_145("Nd-145", "Nd", 25, 85, 60, 145, 8309.1883, "B-", -164.4982, 144.912579151, 1.364);
    static const Species Pm_145("Pm-145", "Pm", 23, 84, 61, 145, 8302.6583, "B-", -616.099, 144.912755748, 3.011);
    static const Species Sm_145("Sm-145", "Sm", 21, 83, 62, 145, 8293.0139, "B-", -2659.8832, 144.913417157, 1.594);
    static const Species Eu_145("Eu-145", "Eu", 19, 82, 63, 145, 8269.2744, "B-", -5064.8984, 144.916272659, 3.285);
    static const Species Gd_145("Gd-145", "Gd", 17, 81, 64, 145, 8228.9485, "B-", -6526.9329, 144.921710051, 21.165);
    static const Species Tb_145("Tb-145", "Tb", 15, 80, 65, 145, 8178.5397, "B-", -8157.0853, 144.928717001, 119.051);
    static const Species Dy_145("Dy-145", "Dy", 13, 79, 66, 145, 8116.8884, "B-", -9122.4943, 144.937473992, 7.0);
    static const Species Ho_145("Ho-145", "Ho", 11, 78, 67, 145, 8048.5792, "B-", -9880.0, 144.947267392, 8.0);
    static const Species Er_145("Er-145", "Er", 9, 77, 68, 145, 7975.0, "B-", -11657.0, 144.957874, 215.0);
    static const Species Tm_145("Tm-145", "Tm", 7, 76, 69, 145, 7889.0, "B-", std::numeric_limits<double>::quiet_NaN(), 144.970389, 210.0);
    static const Species I_146("I-146", "I", 40, 93, 53, 146, 8031.0, "B-", 12415.0, 145.961846, 322.0);
    static const Species Xe_146("Xe-146", "Xe", 38, 92, 54, 146, 8110.4154, "B-", 7355.4298, 145.948518245, 26.0);
    static const Species Cs_146("Cs-146", "Cs", 36, 91, 55, 146, 8155.4365, "B-", 9555.8882, 145.940621867, 3.106);
    static const Species Ba_146("Ba-146", "Ba", 34, 90, 56, 146, 8215.5293, "B-", 4354.9052, 145.9303632, 1.9);
    static const Species La_146("La-146", "La", 32, 89, 57, 146, 8239.9988, "B-", 6404.6947, 145.925688017, 1.797);
    static const Species Ce_146("Ce-146", "Ce", 30, 88, 58, 146, 8278.5081, "B-", 1047.6178, 145.918812294, 15.743);
    static const Species Pr_146("Pr-146", "Pr", 28, 87, 59, 146, 8280.325, "B-", 4252.43, 145.91768763, 36.882);
    static const Species Nd_146("Nd-146", "Nd", 26, 86, 60, 146, 8304.0927, "B-", -1471.556, 145.913122459, 1.366);
    static const Species Pm_146("Pm-146", "Pm", 24, 85, 61, 146, 8288.655, "B-", 1542.0, 145.91470224, 4.589);
    static const Species Sm_146("Sm-146", "Sm", 22, 84, 62, 146, 8293.8581, "B-", -3878.7573, 145.913046835, 3.269);
    static const Species Eu_146("Eu-146", "Eu", 20, 83, 63, 146, 8261.9327, "B-", -1031.7798, 145.917210852, 6.451);
    static const Species Gd_146("Gd-146", "Gd", 18, 82, 64, 146, 8249.5072, "B-", -8322.1791, 145.918318513, 4.376);
    static const Species Tb_146("Tb-146", "Tb", 16, 81, 65, 146, 8187.1474, "B-", -5208.7165, 145.927252739, 48.159);
    static const Species Dy_146("Dy-146", "Dy", 14, 80, 66, 146, 8146.1128, "B-", -11316.7007, 145.932844526, 7.187);
    static const Species Ho_146("Ho-146", "Ho", 12, 79, 67, 146, 8063.2426, "B-", -6916.2074, 145.944993503, 7.071);
    static const Species Er_146("Er-146", "Er", 10, 78, 68, 146, 8010.5127, "B-", -13267.0, 145.952418357, 7.197);
    static const Species Tm_146("Tm-146", "Tm", 8, 77, 69, 146, 7914.0, "B-", std::numeric_limits<double>::quiet_NaN(), 145.966661, 215.0);
    static const Species I_147("I-147", "I", 41, 94, 53, 147, 8001.0, "B-", 11199.0, 146.966505, 322.0);
    static const Species Xe_147("Xe-147", "Xe", 39, 93, 54, 147, 8072.0, "B-", 9520.0, 146.954482, 215.0);
    static const Species Cs_147("Cs-147", "Cs", 37, 92, 55, 147, 8131.801, "B-", 8343.9631, 146.944261512, 9.0);
    static const Species Ba_147("Ba-147", "Ba", 35, 91, 56, 147, 8183.2405, "B-", 6414.3615, 146.9353039, 21.2);
    static const Species La_147("La-147", "La", 33, 90, 57, 147, 8221.5536, "B-", 5335.5046, 146.9284178, 11.5);
    static const Species Ce_147("Ce-147", "Ce", 31, 89, 58, 147, 8252.5274, "B-", 3430.1913, 146.9226899, 9.211);
    static const Species Pr_147("Pr-147", "Pr", 29, 88, 59, 147, 8270.54, "B-", 2702.7013, 146.919007438, 17.02);
    static const Species Nd_147("Nd-147", "Nd", 27, 87, 60, 147, 8283.6036, "B-", 895.1896, 146.916105969, 1.368);
    static const Species Pm_147("Pm-147", "Pm", 25, 86, 61, 147, 8284.3712, "B-", 224.0638, 146.915144944, 1.382);
    static const Species Sm_147("Sm-147", "Sm", 23, 85, 62, 147, 8280.5734, "B-", -1721.4367, 146.914904401, 1.354);
    static const Species Eu_147("Eu-147", "Eu", 21, 84, 63, 147, 8263.5409, "B-", -2187.6833, 146.91675244, 2.758);
    static const Species Gd_147("Gd-147", "Gd", 19, 83, 64, 147, 8243.3366, "B-", -4614.2543, 146.919101014, 2.025);
    static const Species Tb_147("Tb-147", "Tb", 17, 82, 65, 147, 8206.625, "B-", -6546.6266, 146.92405462, 8.691);
    static const Species Dy_147("Dy-147", "Dy", 15, 81, 66, 147, 8156.768, "B-", -8438.946, 146.931082712, 9.5);
    static const Species Ho_147("Ho-147", "Ho", 13, 80, 67, 147, 8094.0381, "B-", -9149.2869, 146.940142293, 5.368);
    static const Species Er_147("Er-147", "Er", 11, 79, 68, 147, 8026.476, "B-", -10633.4071, 146.949964456, 41.0);
    static const Species Tm_147("Tm-147", "Tm", 9, 78, 69, 147, 7948.8178, "B-", std::numeric_limits<double>::quiet_NaN(), 146.961379887, 7.341);
    static const Species Xe_148("Xe-148", "Xe", 40, 94, 54, 148, 8047.0, "B-", 8261.0, 147.958508, 322.0);
    static const Species Cs_148("Cs-148", "Cs", 38, 93, 55, 148, 8097.5469, "B-", 10633.9614, 147.949639026, 14.0);
    static const Species Ba_148("Ba-148", "Ba", 36, 92, 56, 148, 8164.1118, "B-", 5163.8307, 147.938223, 1.6);
    static const Species La_148("La-148", "La", 34, 91, 57, 148, 8193.7165, "B-", 7689.6824, 147.9326794, 20.9);
    static const Species Ce_148("Ce-148", "Ce", 32, 90, 58, 148, 8240.3876, "B-", 2137.0282, 147.924424186, 12.017);
    static const Species Pr_148("Pr-148", "Pr", 30, 89, 59, 148, 8249.5409, "B-", 4872.6133, 147.922129992, 16.147);
    static const Species Nd_148("Nd-148", "Nd", 28, 88, 60, 148, 8277.1778, "B-", -542.1891, 147.916899027, 2.203);
    static const Species Pm_148("Pm-148", "Pm", 26, 87, 61, 148, 8268.2283, "B-", 2470.1898, 147.917481091, 6.108);
    static const Species Sm_148("Sm-148", "Sm", 24, 86, 62, 148, 8279.6326, "B-", -3038.5844, 147.914829233, 1.337);
    static const Species Eu_148("Eu-148", "Eu", 22, 85, 63, 148, 8253.8155, "B-", -28.063, 147.918091288, 10.693);
    static const Species Gd_148("Gd-148", "Gd", 20, 84, 64, 148, 8248.3398, "B-", -5732.4723, 147.918121414, 1.566);
    static const Species Tb_148("Tb-148", "Tb", 18, 83, 65, 148, 8204.3207, "B-", -2677.5501, 147.924275476, 13.379);
    static const Species Dy_148("Dy-148", "Dy", 16, 82, 66, 148, 8180.943, "B-", -9868.2304, 147.927149944, 9.365);
    static const Species Ho_148("Ho-148", "Ho", 14, 81, 67, 148, 8108.9797, "B-", -6512.1694, 147.937743925, 90.0);
    static const Species Er_148("Er-148", "Er", 12, 80, 68, 148, 8059.6924, "B-", -12713.963, 147.944735026, 11.0);
    static const Species Tm_148("Tm-148", "Tm", 10, 79, 69, 148, 7968.5011, "B-", -8535.0, 147.958384026, 11.0);
    static const Species Yb_148("Yb-148", "Yb", 8, 78, 70, 148, 7906.0, "B-", std::numeric_limits<double>::quiet_NaN(), 147.967547, 429.0);
    static const Species Xe_149("Xe-149", "Xe", 41, 95, 54, 149, 8009.0, "B-", 10300.0, 148.964573, 322.0);
    static const Species Cs_149("Cs-149", "Cs", 39, 94, 55, 149, 8073.0, "B-", 9531.0, 148.953516, 429.0);
    static const Species Ba_149("Ba-149", "Ba", 37, 93, 56, 149, 8131.8495, "B-", 7389.301, 148.943284, 2.7);
    static const Species La_149("La-149", "La", 35, 92, 57, 149, 8176.1915, "B-", 6450.0, 148.935351259, 214.99);
    static const Species Ce_149("Ce-149", "Ce", 33, 91, 58, 149, 8214.2294, "B-", 4369.4525, 148.9284269, 11.0);
    static const Species Pr_149("Pr-149", "Pr", 31, 90, 59, 149, 8238.304, "B-", 3336.1617, 148.9237361, 10.6);
    static const Species Nd_149("Nd-149", "Nd", 29, 89, 60, 149, 8255.4437, "B-", 1688.8697, 148.920154583, 2.205);
    static const Species Pm_149("Pm-149", "Pm", 27, 88, 61, 149, 8261.5277, "B-", 1071.4936, 148.918341507, 2.344);
    static const Species Sm_149("Sm-149", "Sm", 25, 87, 62, 149, 8263.4683, "B-", -694.5812, 148.917191211, 1.241);
    static const Species Eu_149("Eu-149", "Eu", 23, 86, 63, 149, 8253.556, "B-", -1314.1437, 148.917936875, 4.19);
    static const Species Gd_149("Gd-149", "Gd", 21, 85, 64, 149, 8239.4856, "B-", -3638.5336, 148.919347666, 3.553);
    static const Species Tb_149("Tb-149", "Tb", 19, 84, 65, 149, 8209.8153, "B-", -3794.6493, 148.923253792, 3.895);
    static const Species Dy_149("Dy-149", "Dy", 17, 83, 66, 149, 8179.0972, "B-", -6048.1366, 148.927327516, 9.858);
    static const Species Ho_149("Ho-149", "Ho", 15, 82, 67, 149, 8133.255, "B-", -7904.2328, 148.933820457, 12.866);
    static const Species Er_149("Er-149", "Er", 13, 81, 68, 149, 8074.9558, "B-", -9801.0, 148.942306, 30.0);
    static const Species Tm_149("Tm-149", "Tm", 11, 80, 69, 149, 8004.0, "B-", -10611.0, 148.952828, 215.0);
    static const Species Yb_149("Yb-149", "Yb", 9, 79, 70, 149, 7927.0, "B-", std::numeric_limits<double>::quiet_NaN(), 148.964219, 322.0);
    static const Species Xe_150("Xe-150", "Xe", 42, 96, 54, 150, 7983.0, "B-", 9180.0, 149.968878, 322.0);
    static const Species Cs_150("Cs-150", "Cs", 40, 95, 55, 150, 8039.0, "B-", 11720.0, 149.959023, 429.0);
    static const Species Ba_150("Ba-150", "Ba", 38, 94, 56, 150, 8111.8405, "B-", 6421.3477, 149.9464411, 6.1);
    static const Species La_150("La-150", "La", 36, 93, 57, 150, 8149.4339, "B-", 8535.7155, 149.9395475, 2.7);
    static const Species Ce_150("Ce-150", "Ce", 34, 92, 58, 150, 8201.123, "B-", 3453.6462, 149.930384032, 12.556);
    static const Species Pr_150("Pr-150", "Pr", 32, 91, 59, 150, 8218.9316, "B-", 5379.4422, 149.926676391, 9.677);
    static const Species Nd_150("Nd-150", "Nd", 30, 90, 60, 150, 8249.5789, "B-", -82.6155, 149.920901322, 1.211);
    static const Species Pm_150("Pm-150", "Pm", 28, 89, 61, 150, 8243.8125, "B-", 3454.0, 149.920990014, 21.504);
    static const Species Sm_150("Sm-150", "Sm", 26, 88, 62, 150, 8261.6235, "B-", -2258.9662, 149.917281993, 1.193);
    static const Species Eu_150("Eu-150", "Eu", 24, 87, 63, 150, 8241.3481, "B-", 971.6815, 149.919707092, 6.688);
    static const Species Gd_150("Gd-150", "Gd", 22, 86, 64, 150, 8242.6103, "B-", -4658.2621, 149.918663949, 6.5);
    static const Species Tb_150("Tb-150", "Tb", 20, 85, 65, 150, 8206.3396, "B-", -1796.1707, 149.923664799, 7.912);
    static const Species Dy_150("Dy-150", "Dy", 18, 84, 66, 150, 8189.1495, "B-", -7363.7264, 149.925593068, 4.636);
    static const Species Ho_150("Ho-150", "Ho", 16, 83, 67, 150, 8134.8423, "B-", -4114.5689, 149.933498353, 15.209);
    static const Species Er_150("Er-150", "Er", 14, 82, 68, 150, 8102.1962, "B-", -11340.0, 149.937915524, 18.458);
    static const Species Tm_150("Tm-150", "Tm", 12, 81, 69, 150, 8021.0, "B-", -7661.0, 149.95009, 210.0);
    static const Species Yb_150("Yb-150", "Yb", 10, 80, 70, 150, 7965.0, "B-", -14059.0, 149.958314, 322.0);
    static const Species Lu_150("Lu-150", "Lu", 8, 79, 71, 150, 7866.0, "B-", std::numeric_limits<double>::quiet_NaN(), 149.973407, 322.0);
    static const Species Cs_151("Cs-151", "Cs", 41, 96, 55, 151, 8013.0, "B-", 10660.0, 150.963199, 537.0);
    static const Species Ba_151("Ba-151", "Ba", 39, 95, 56, 151, 8079.0, "B-", 8370.0, 150.951755, 429.0);
    static const Species La_151("La-151", "La", 37, 94, 57, 151, 8129.0436, "B-", 7914.7191, 150.942769, 467.5);
    static const Species Ce_151("Ce-151", "Ce", 35, 93, 58, 151, 8176.2779, "B-", 5554.6233, 150.9342722, 19.0);
    static const Species Pr_151("Pr-151", "Pr", 33, 92, 59, 151, 8207.8824, "B-", 4163.5021, 150.928309066, 12.506);
    static const Species Nd_151("Nd-151", "Nd", 31, 91, 60, 151, 8230.2741, "B-", 2443.0767, 150.923839363, 1.215);
    static const Species Pm_151("Pm-151", "Pm", 29, 90, 61, 151, 8241.2723, "B-", 1190.2198, 150.921216613, 4.949);
    static const Species Sm_151("Sm-151", "Sm", 27, 89, 62, 151, 8243.9735, "B-", 76.6182, 150.919938859, 1.191);
    static const Species Eu_151("Eu-151", "Eu", 25, 88, 63, 151, 8239.2998, "B-", -464.1779, 150.919856606, 1.251);
    static const Species Gd_151("Gd-151", "Gd", 23, 87, 64, 151, 8231.0446, "B-", -2565.3796, 150.920354922, 3.212);
    static const Species Tb_151("Tb-151", "Tb", 21, 86, 65, 151, 8208.8743, "B-", -2871.1525, 150.92310897, 4.395);
    static const Species Dy_151("Dy-151", "Dy", 19, 85, 66, 151, 8184.6789, "B-", -5129.6421, 150.926191279, 3.486);
    static const Species Ho_151("Ho-151", "Ho", 17, 84, 67, 151, 8145.5267, "B-", -5356.4558, 150.931698176, 8.908);
    static const Species Er_151("Er-151", "Er", 15, 83, 68, 151, 8104.8723, "B-", -7494.6768, 150.937448567, 17.681);
    static const Species Tm_151("Tm-151", "Tm", 13, 82, 69, 151, 8050.0576, "B-", -9229.2615, 150.945494433, 20.799);
    static const Species Yb_151("Yb-151", "Yb", 11, 81, 70, 151, 7983.7556, "B-", -11242.0, 150.955402453, 322.591);
    static const Species Lu_151("Lu-151", "Lu", 9, 80, 71, 151, 7904.0, "B-", std::numeric_limits<double>::quiet_NaN(), 150.967471, 322.0);
    static const Species Cs_152("Cs-152", "Cs", 42, 97, 55, 152, 7980.0, "B-", 12480.0, 151.968728, 537.0);
    static const Species Ba_152("Ba-152", "Ba", 40, 96, 56, 152, 8057.0, "B-", 7680.0, 151.95533, 429.0);
    static const Species La_152("La-152", "La", 38, 95, 57, 152, 8102.0, "B-", 9690.0, 151.947085, 322.0);
    static const Species Ce_152("Ce-152", "Ce", 36, 94, 58, 152, 8161.0, "B-", 4778.0, 151.936682, 215.0);
    static const Species Pr_152("Pr-152", "Pr", 34, 93, 59, 152, 8187.1049, "B-", 6391.5934, 151.9315529, 19.9);
    static const Species Nd_152("Nd-152", "Nd", 32, 92, 60, 152, 8224.0078, "B-", 1104.805, 151.924691242, 26.276);
    static const Species Pm_152("Pm-152", "Pm", 30, 91, 61, 152, 8226.1293, "B-", 3508.5089, 151.923505185, 27.809);
    static const Species Sm_152("Sm-152", "Sm", 28, 90, 62, 152, 8244.0645, "B-", -1874.4774, 151.919738646, 1.09);
    static const Species Eu_152("Eu-152", "Eu", 26, 89, 63, 152, 8226.5854, "B-", 1818.8037, 151.92175098, 1.252);
    static const Species Gd_152("Gd-152", "Gd", 24, 88, 64, 152, 8233.4042, "B-", -3990.0, 151.919798414, 1.081);
    static const Species Tb_152("Tb-152", "Tb", 22, 87, 65, 152, 8202.0072, "B-", -599.3405, 151.924081855, 42.955);
    static const Species Dy_152("Dy-152", "Dy", 20, 86, 66, 152, 8192.9171, "B-", -6513.3275, 151.924725274, 4.93);
    static const Species Ho_152("Ho-152", "Ho", 18, 85, 67, 152, 8144.9193, "B-", -3104.4174, 151.931717618, 13.449);
    static const Species Er_152("Er-152", "Er", 16, 84, 68, 152, 8119.3485, "B-", -8779.9397, 151.935050347, 9.478);
    static const Species Tm_152("Tm-152", "Tm", 14, 83, 69, 152, 8056.4387, "B-", -5449.8923, 151.944476, 58.0);
    static const Species Yb_152("Yb-152", "Yb", 12, 82, 70, 152, 8015.4371, "B-", -12848.0, 151.950326699, 160.718);
    static const Species Lu_152("Lu-152", "Lu", 10, 81, 71, 152, 7926.0, "B-", std::numeric_limits<double>::quiet_NaN(), 151.96412, 210.0);
    static const Species Ba_153("Ba-153", "Ba", 41, 97, 56, 153, 8023.0, "B-", 9590.0, 152.960848, 429.0);
    static const Species La_153("La-153", "La", 39, 96, 57, 153, 8081.0, "B-", 8850.0, 152.950553, 322.0);
    static const Species Ce_153("Ce-153", "Ce", 37, 95, 58, 153, 8134.0, "B-", 6659.0, 152.941052, 215.0);
    static const Species Pr_153("Pr-153", "Pr", 35, 94, 59, 153, 8172.0371, "B-", 5761.8901, 152.933903511, 12.755);
    static const Species Nd_153("Nd-153", "Nd", 33, 93, 60, 153, 8204.5832, "B-", 3317.6236, 152.927717868, 2.949);
    static const Species Pm_153("Pm-153", "Pm", 31, 92, 61, 153, 8221.1536, "B-", 1912.0559, 152.924156252, 9.729);
    static const Species Sm_153("Sm-153", "Sm", 29, 91, 62, 153, 8228.5373, "B-", 807.4073, 152.922103576, 1.1);
    static const Species Eu_153("Eu-153", "Eu", 27, 90, 63, 153, 8228.7011, "B-", -484.5225, 152.921236789, 1.257);
    static const Species Gd_153("Gd-153", "Gd", 25, 89, 64, 153, 8220.4209, "B-", -1569.334, 152.921756945, 1.075);
    static const Species Tb_153("Tb-153", "Tb", 23, 88, 65, 153, 8205.0504, "B-", -2170.4134, 152.923441694, 4.237);
    static const Species Dy_153("Dy-153", "Dy", 21, 87, 66, 153, 8185.7514, "B-", -4131.1229, 152.925771729, 4.295);
    static const Species Ho_153("Ho-153", "Ho", 19, 86, 67, 153, 8153.6372, "B-", -4545.3918, 152.930206671, 5.438);
    static const Species Er_153("Er-153", "Er", 17, 85, 68, 153, 8118.8154, "B-", -6494.0723, 152.93508635, 9.967);
    static const Species Tm_153("Tm-153", "Tm", 15, 84, 69, 153, 8071.2571, "B-", -6813.0, 152.942058023, 12.86);
    static const Species Yb_153("Yb-153", "Yb", 13, 83, 70, 153, 8022.0, "B-", -8784.0, 152.949372, 215.0);
    static const Species Lu_153("Lu-153", "Lu", 11, 82, 71, 153, 7959.0882, "B-", -11075.0, 152.958802248, 161.05);
    static const Species Hf_153("Hf-153", "Hf", 9, 81, 72, 153, 7882.0, "B-", std::numeric_limits<double>::quiet_NaN(), 152.970692, 322.0);
    static const Species Ba_154("Ba-154", "Ba", 42, 98, 56, 154, 8001.0, "B-", 8610.0, 153.964659, 537.0);
    static const Species La_154("La-154", "La", 40, 97, 57, 154, 8051.0, "B-", 10690.0, 153.955416, 322.0);
    static const Species Ce_154("Ce-154", "Ce", 38, 96, 58, 154, 8116.0, "B-", 5640.0, 153.94394, 215.0);
    static const Species Pr_154("Pr-154", "Pr", 36, 95, 59, 154, 8147.2994, "B-", 7720.0, 153.937885165, 107.36);
    static const Species Nd_154("Nd-154", "Nd", 34, 94, 60, 154, 8192.3491, "B-", 2687.0, 153.929597404, 1.1);
    static const Species Pm_154("Pm-154", "Pm", 32, 93, 61, 154, 8204.717, "B-", 4188.9614, 153.926712791, 26.861);
    static const Species Sm_154("Sm-154", "Sm", 30, 92, 62, 154, 8226.8379, "B-", -717.1969, 153.922215756, 1.4);
    static const Species Eu_154("Eu-154", "Eu", 28, 91, 63, 154, 8217.1006, "B-", 1967.9913, 153.922985699, 1.275);
    static const Species Gd_154("Gd-154", "Gd", 26, 90, 64, 154, 8224.7996, "B-", -3549.6514, 153.920872974, 1.066);
    static const Species Tb_154("Tb-154", "Tb", 24, 89, 65, 154, 8196.6697, "B-", 237.308, 153.924683681, 48.641);
    static const Species Dy_154("Dy-154", "Dy", 22, 88, 66, 154, 8193.1305, "B-", -5754.6363, 153.92442892, 7.977);
    static const Species Ho_154("Ho-154", "Ho", 20, 87, 67, 154, 8150.6825, "B-", -2034.4045, 153.930606776, 8.82);
    static const Species Er_154("Er-154", "Er", 18, 86, 68, 154, 8132.3919, "B-", -8177.8316, 153.932790799, 5.325);
    static const Species Tm_154("Tm-154", "Tm", 16, 85, 69, 154, 8074.209, "B-", -4495.0495, 153.941570062, 15.471);
    static const Species Yb_154("Yb-154", "Yb", 14, 84, 70, 154, 8039.9402, "B-", -10265.0, 153.946395696, 18.551);
    static const Species Lu_154("Lu-154", "Lu", 12, 83, 71, 154, 7968.0, "B-", -6937.0, 153.957416, 216.0);
    static const Species Hf_154("Hf-154", "Hf", 10, 82, 72, 154, 7918.0, "B-", std::numeric_limits<double>::quiet_NaN(), 153.964863, 322.0);
    static const Species La_155("La-155", "La", 41, 98, 57, 155, 8028.0, "B-", 9850.0, 154.95928, 429.0);
    static const Species Ce_155("Ce-155", "Ce", 39, 97, 58, 155, 8087.0, "B-", 7635.0, 154.948706, 322.0);
    static const Species Pr_155("Pr-155", "Pr", 37, 96, 59, 155, 8131.0398, "B-", 6868.4609, 154.940509193, 18.462);
    static const Species Nd_155("Nd-155", "Nd", 35, 95, 60, 155, 8170.305, "B-", 4656.2095, 154.933135598, 9.826);
    static const Species Pm_155("Pm-155", "Pm", 33, 94, 61, 155, 8195.2977, "B-", 3251.1999, 154.928136951, 5.065);
    static const Species Sm_155("Sm-155", "Sm", 31, 93, 62, 155, 8211.2258, "B-", 1627.1314, 154.924646645, 1.429);
    static const Species Eu_155("Eu-155", "Eu", 29, 92, 63, 155, 8216.676, "B-", 251.9612, 154.922899847, 1.343);
    static const Species Gd_155("Gd-155", "Gd", 27, 91, 64, 155, 8213.2541, "B-", -819.8588, 154.922629356, 1.055);
    static const Species Tb_155("Tb-155", "Tb", 25, 90, 65, 155, 8202.9173, "B-", -2094.5, 154.923509511, 10.552);
    static const Species Dy_155("Dy-155", "Dy", 23, 89, 66, 155, 8184.357, "B-", -3116.1405, 154.925758049, 10.354);
    static const Species Ho_155("Ho-155", "Ho", 21, 88, 67, 155, 8159.2055, "B-", -3830.6268, 154.929103363, 18.754);
    static const Species Er_155("Er-155", "Er", 19, 87, 68, 155, 8129.4444, "B-", -5583.2509, 154.93321571, 6.52);
    static const Species Tm_155("Tm-155", "Tm", 17, 86, 69, 155, 8088.376, "B-", -6123.3072, 154.939209576, 10.651);
    static const Species Yb_155("Yb-155", "Yb", 15, 85, 70, 155, 8043.8234, "B-", -7957.5578, 154.945783216, 17.82);
    static const Species Lu_155("Lu-155", "Lu", 13, 84, 71, 155, 7987.4369, "B-", -8235.0, 154.954326005, 20.66);
    static const Species Hf_155("Hf-155", "Hf", 11, 83, 72, 155, 7929.0, "B-", -10322.0, 154.963167, 322.0);
    static const Species Ta_155("Ta-155", "Ta", 9, 82, 73, 155, 7858.0, "B-", std::numeric_limits<double>::quiet_NaN(), 154.974248, 322.0);
    static const Species La_156("La-156", "La", 42, 99, 57, 156, 7997.0, "B-", 11769.0, 155.964519, 429.0);
    static const Species Ce_156("Ce-156", "Ce", 40, 98, 58, 156, 8068.0, "B-", 6630.0, 155.951884, 322.0);
    static const Species Pr_156("Pr-156", "Pr", 38, 97, 59, 156, 8105.2337, "B-", 8752.8229, 155.9447669, 1.1);
    static const Species Nd_156("Nd-156", "Nd", 36, 96, 60, 156, 8156.3265, "B-", 3964.7175, 155.935370358, 1.4);
    static const Species Pm_156("Pm-156", "Pm", 34, 95, 61, 156, 8176.7263, "B-", 5193.8878, 155.931114059, 1.275);
    static const Species Sm_156("Sm-156", "Sm", 32, 94, 62, 156, 8205.0054, "B-", 722.109, 155.925538191, 9.148);
    static const Species Eu_156("Eu-156", "Eu", 30, 93, 63, 156, 8204.6192, "B-", 2452.4891, 155.924762976, 3.791);
    static const Species Gd_156("Gd-156", "Gd", 28, 92, 64, 156, 8215.3253, "B-", -2444.323, 155.92213012, 1.054);
    static const Species Tb_156("Tb-156", "Tb", 26, 91, 65, 156, 8194.6415, "B-", 438.3762, 155.924754209, 4.044);
    static const Species Dy_156("Dy-156", "Dy", 24, 90, 66, 156, 8192.4366, "B-", -4990.9836, 155.924283593, 1.06);
    static const Species Ho_156("Ho-156", "Ho", 22, 89, 67, 156, 8155.428, "B-", -1326.7201, 155.929641634, 41.249);
    static const Species Er_156("Er-156", "Er", 20, 88, 68, 156, 8141.9084, "B-", -7377.2657, 155.931065926, 26.44);
    static const Species Tm_156("Tm-156", "Tm", 18, 87, 69, 156, 8089.6032, "B-", -3568.8794, 155.938985746, 15.328);
    static const Species Yb_156("Yb-156", "Yb", 16, 86, 70, 156, 8061.7107, "B-", -9565.988, 155.942817096, 9.992);
    static const Species Lu_156("Lu-156", "Lu", 14, 85, 71, 156, 7995.3752, "B-", -5880.0352, 155.953086606, 58.102);
    static const Species Hf_156("Hf-156", "Hf", 12, 84, 72, 156, 7952.6676, "B-", -11819.0, 155.959399083, 160.752);
    static const Species Ta_156("Ta-156", "Ta", 10, 83, 73, 156, 7872.0, "B-", std::numeric_limits<double>::quiet_NaN(), 155.972087, 322.0);
    static const Species La_157("La-157", "La", 43, 100, 57, 157, 7972.0, "B-", 10860.0, 156.968792, 322.0);
    static const Species Ce_157("Ce-157", "Ce", 41, 99, 58, 157, 8037.0, "B-", 8504.0, 156.957133, 429.0);
    static const Species Pr_157("Pr-157", "Pr", 39, 98, 59, 157, 8085.817, "B-", 8059.3109, 156.9480031, 3.4);
    static const Species Nd_157("Nd-157", "Nd", 37, 97, 60, 157, 8132.1671, "B-", 5802.9994, 156.939351074, 2.294);
    static const Species Pm_157("Pm-157", "Pm", 35, 96, 61, 157, 8164.1458, "B-", 4380.5376, 156.933121298, 7.521);
    static const Species Sm_157("Sm-157", "Sm", 33, 95, 62, 157, 8187.0642, "B-", 2781.4807, 156.928418598, 4.759);
    static const Species Eu_157("Eu-157", "Eu", 31, 94, 63, 157, 8199.7975, "B-", 1364.7614, 156.925432556, 4.545);
    static const Species Gd_157("Gd-157", "Gd", 29, 93, 64, 157, 8203.5072, "B-", -60.0473, 156.923967424, 1.048);
    static const Species Tb_157("Tb-157", "Tb", 27, 92, 65, 157, 8198.1416, "B-", -1339.1791, 156.924031888, 1.092);
    static const Species Dy_157("Dy-157", "Dy", 25, 91, 66, 157, 8184.6287, "B-", -2591.8068, 156.925469555, 5.532);
    static const Species Ho_157("Ho-157", "Ho", 23, 90, 67, 157, 8163.1373, "B-", -3419.2146, 156.928251974, 25.194);
    static const Species Er_157("Er-157", "Er", 21, 89, 68, 157, 8136.3757, "B-", -4704.369, 156.931922652, 28.454);
    static const Species Tm_157("Tm-157", "Tm", 19, 88, 69, 157, 8101.4285, "B-", -5289.3667, 156.936973, 30.0);
    static const Species Yb_157("Yb-157", "Yb", 17, 87, 70, 157, 8062.7551, "B-", -6980.0942, 156.942651368, 11.706);
    static const Species Lu_157("Lu-157", "Lu", 15, 86, 71, 157, 8013.3128, "B-", -7585.0, 156.950144807, 12.961);
    static const Species Hf_157("Hf-157", "Hf", 13, 85, 72, 157, 7960.0, "B-", -9259.0, 156.958288, 215.0);
    static const Species Ta_157("Ta-157", "Ta", 11, 84, 73, 157, 7896.0608, "B-", -9906.0, 156.968227445, 161.087);
    static const Species W_157("W-157", "W", 9, 83, 74, 157, 7828.0, "B-", std::numeric_limits<double>::quiet_NaN(), 156.978862, 429.0);
    static const Species Ce_158("Ce-158", "Ce", 42, 100, 58, 158, 8015.0, "B-", 7610.0, 157.960773, 429.0);
    static const Species Pr_158("Pr-158", "Pr", 40, 99, 59, 158, 8059.0, "B-", 9685.0, 157.952603, 322.0);
    static const Species Nd_158("Nd-158", "Nd", 38, 98, 60, 158, 8114.9529, "B-", 5271.0196, 157.94220562, 1.4);
    static const Species Pm_158("Pm-158", "Pm", 36, 97, 61, 158, 8143.3622, "B-", 6145.7062, 157.936546948, 0.953);
    static const Species Sm_158("Sm-158", "Sm", 34, 96, 62, 158, 8177.3075, "B-", 2018.6123, 157.929949262, 5.133);
    static const Species Eu_158("Eu-158", "Eu", 32, 95, 63, 158, 8185.132, "B-", 3419.5081, 157.927782192, 2.181);
    static const Species Gd_158("Gd-158", "Gd", 30, 94, 64, 158, 8201.8229, "B-", -1219.0862, 157.9241112, 1.048);
    static const Species Tb_158("Tb-158", "Tb", 28, 93, 65, 158, 8189.1556, "B-", 936.2686, 157.925419942, 1.36);
    static const Species Dy_158("Dy-158", "Dy", 26, 92, 66, 158, 8190.1298, "B-", -4219.7555, 157.924414817, 2.509);
    static const Species Ho_158("Ho-158", "Ho", 24, 91, 67, 158, 8158.4709, "B-", -883.5812, 157.92894491, 29.099);
    static const Species Er_158("Er-158", "Er", 22, 90, 68, 158, 8147.927, "B-", -6600.6151, 157.929893474, 27.074);
    static const Species Tm_158("Tm-158", "Tm", 20, 89, 69, 158, 8101.1994, "B-", -2693.5796, 157.936979525, 27.074);
    static const Species Yb_158("Yb-158", "Yb", 18, 88, 70, 158, 8079.1999, "B-", -8797.4201, 157.939871202, 8.559);
    static const Species Lu_158("Lu-158", "Lu", 16, 87, 71, 158, 8018.5685, "B-", -5109.801, 157.94931562, 16.236);
    static const Species Hf_158("Hf-158", "Hf", 14, 86, 72, 158, 7981.2764, "B-", -10984.0, 157.954801217, 18.78);
    static const Species Ta_158("Ta-158", "Ta", 12, 85, 73, 158, 7907.0, "B-", -7426.0, 157.966593, 215.0);
    static const Species W_158("W-158", "W", 10, 84, 74, 158, 7855.0, "B-", std::numeric_limits<double>::quiet_NaN(), 157.974565, 322.0);
    static const Species Ce_159("Ce-159", "Ce", 43, 101, 58, 159, 7983.0, "B-", 9430.0, 158.966355, 537.0);
    static const Species Pr_159("Pr-159", "Pr", 41, 100, 59, 159, 8037.0, "B-", 8954.0, 158.956232, 429.0);
    static const Species Nd_159("Nd-159", "Nd", 39, 99, 60, 159, 8088.8224, "B-", 6830.3441, 158.946619085, 32.0);
    static const Species Pm_159("Pm-159", "Pm", 37, 98, 61, 159, 8126.8601, "B-", 5653.4982, 158.939286409, 10.777);
    static const Species Sm_159("Sm-159", "Sm", 35, 97, 62, 159, 8157.4963, "B-", 3835.5363, 158.93321713, 6.37);
    static const Species Eu_159("Eu-159", "Eu", 33, 96, 63, 159, 8176.6987, "B-", 2518.4717, 158.929099512, 4.637);
    static const Species Gd_159("Gd-159", "Gd", 31, 95, 64, 159, 8187.6177, "B-", 970.7242, 158.926395822, 1.051);
    static const Species Tb_159("Tb-159", "Tb", 29, 94, 65, 159, 8188.8025, "B-", -365.3613, 158.925353707, 1.184);
    static const Species Dy_159("Dy-159", "Dy", 27, 93, 66, 159, 8181.5842, "B-", -1837.6, 158.925745938, 1.544);
    static const Species Ho_159("Ho-159", "Ho", 25, 92, 67, 159, 8165.1066, "B-", -2768.5, 158.927718683, 3.268);
    static const Species Er_159("Er-159", "Er", 23, 91, 68, 159, 8142.7742, "B-", -3990.7162, 158.93069079, 3.91);
    static const Species Tm_159("Tm-159", "Tm", 21, 90, 69, 159, 8112.7549, "B-", -4736.8869, 158.934975, 30.0);
    static const Species Yb_159("Yb-159", "Yb", 19, 89, 70, 159, 8078.0428, "B-", -6124.9081, 158.940060257, 18.874);
    static const Species Lu_159("Lu-159", "Lu", 17, 88, 71, 159, 8034.6009, "B-", -6856.003, 158.946635615, 40.433);
    static const Species Hf_159("Hf-159", "Hf", 15, 87, 72, 159, 7986.561, "B-", -8413.4492, 158.953995837, 18.049);
    static const Species Ta_159("Ta-159", "Ta", 13, 86, 73, 159, 7928.7258, "B-", -9005.0, 158.963028046, 21.137);
    static const Species W_159("W-159", "W", 11, 85, 74, 159, 7867.0, "B-", -10629.0, 158.972696, 322.0);
    static const Species Re_159("Re-159", "Re", 9, 84, 75, 159, 7795.0, "B-", std::numeric_limits<double>::quiet_NaN(), 158.984106, 327.0);
    static const Species Pr_160("Pr-160", "Pr", 42, 101, 59, 160, 8009.0, "B-", 10525.0, 159.961138, 429.0);
    static const Species Nd_160("Nd-160", "Nd", 40, 100, 60, 160, 8069.9662, "B-", 6170.1238, 159.949839172, 50.0);
    static const Species Pm_160("Pm-160", "Pm", 38, 99, 61, 160, 8103.6398, "B-", 7338.5338, 159.943215272, 2.2);
    static const Species Sm_160("Sm-160", "Sm", 36, 98, 62, 160, 8144.6159, "B-", 3260.2763, 159.935337032, 2.1);
    static const Species Eu_160("Eu-160", "Eu", 34, 97, 63, 160, 8160.103, "B-", 4448.6112, 159.931836982, 0.97);
    static const Species Gd_160("Gd-160", "Gd", 32, 96, 64, 160, 8183.0171, "B-", -105.5863, 159.927061202, 1.206);
    static const Species Tb_160("Tb-160", "Tb", 30, 95, 65, 160, 8177.4676, "B-", 1835.9516, 159.927174553, 1.191);
    static const Species Dy_160("Dy-160", "Dy", 28, 94, 66, 160, 8184.0526, "B-", -3290.0, 159.925203578, 0.751);
    static const Species Ho_160("Ho-160", "Ho", 26, 93, 67, 160, 8158.6004, "B-", -318.2488, 159.928735538, 16.12);
    static const Species Er_160("Er-160", "Er", 24, 92, 68, 160, 8151.7217, "B-", -5763.1395, 159.929077193, 26.029);
    static const Species Tm_160("Tm-160", "Tm", 22, 91, 69, 160, 8110.8124, "B-", -2137.8101, 159.935264177, 35.089);
    static const Species Yb_160("Yb-160", "Yb", 20, 90, 70, 160, 8092.5614, "B-", -7893.2846, 159.93755921, 5.9);
    static const Species Lu_160("Lu-160", "Lu", 18, 89, 71, 160, 8038.3387, "B-", -4331.1951, 159.946033, 61.0);
    static const Species Hf_160("Hf-160", "Hf", 16, 88, 72, 160, 8006.3791, "B-", -10115.0475, 159.950682728, 10.241);
    static const Species Ta_160("Ta-160", "Ta", 14, 87, 73, 160, 7938.2704, "B-", -6494.6267, 159.961541678, 58.31);
    static const Species W_160("W-160", "W", 12, 86, 74, 160, 7892.7893, "B-", -12451.0, 159.968513946, 160.828);
    static const Species Re_160("Re-160", "Re", 10, 85, 75, 160, 7810.0, "B-", std::numeric_limits<double>::quiet_NaN(), 159.98188, 322.0);
    static const Species Pr_161("Pr-161", "Pr", 43, 102, 59, 161, 7986.0, "B-", 9741.0, 160.965121, 537.0);
    static const Species Nd_161("Nd-161", "Nd", 41, 101, 60, 161, 8042.0, "B-", 7856.0, 160.954664, 429.0);
    static const Species Pm_161("Pm-161", "Pm", 39, 100, 61, 161, 8085.9977, "B-", 6585.4546, 160.946229837, 9.7);
    static const Species Sm_161("Sm-161", "Sm", 37, 99, 62, 161, 8122.0418, "B-", 5119.5577, 160.939160062, 7.318);
    static const Species Eu_161("Eu-161", "Eu", 35, 98, 63, 161, 8148.981, "B-", 3714.5406, 160.933663991, 11.164);
    static const Species Gd_161("Gd-161", "Gd", 33, 97, 64, 161, 8167.1934, "B-", 1955.6357, 160.929676267, 1.614);
    static const Species Tb_161("Tb-161", "Tb", 31, 96, 65, 161, 8174.4809, "B-", 593.7166, 160.927576806, 1.308);
    static const Species Dy_161("Dy-161", "Dy", 29, 95, 66, 161, 8173.3093, "B-", -859.2003, 160.926939425, 0.748);
    static const Species Ho_161("Ho-161", "Ho", 27, 94, 67, 161, 8163.1134, "B-", -1994.9954, 160.927861815, 2.309);
    static const Species Er_161("Er-161", "Er", 25, 93, 68, 161, 8145.8628, "B-", -3302.5839, 160.93000353, 9.419);
    static const Species Tm_161("Tm-161", "Tm", 23, 92, 69, 161, 8120.4906, "B-", -4064.4665, 160.933549, 30.0);
    static const Species Yb_161("Yb-161", "Yb", 21, 91, 70, 161, 8090.3861, "B-", -5271.8989, 160.937912384, 16.211);
    static const Species Lu_161("Lu-161", "Lu", 19, 90, 71, 161, 8052.7821, "B-", -6246.5319, 160.943572, 30.0);
    static const Species Hf_161("Hf-161", "Hf", 17, 89, 72, 161, 8009.1245, "B-", -7537.2421, 160.950277927, 25.174);
    static const Species Ta_161("Ta-161", "Ta", 15, 88, 73, 161, 7957.45, "B-", -8272.0, 160.958369489, 26.174);
    static const Species W_161("W-161", "W", 13, 87, 74, 161, 7901.0, "B-", -9664.0, 160.967249, 215.0);
    static const Species Re_161("Re-161", "Re", 11, 86, 75, 161, 7836.3292, "B-", -10647.0, 160.977624313, 160.93);
    static const Species Os_161("Os-161", "Os", 9, 85, 76, 161, 7765.0, "B-", std::numeric_limits<double>::quiet_NaN(), 160.989054, 429.0);
    static const Species Nd_162("Nd-162", "Nd", 42, 102, 60, 162, 8022.0, "B-", 7030.0, 161.958121, 429.0);
    static const Species Pm_162("Pm-162", "Pm", 40, 101, 61, 162, 8061.0, "B-", 8339.0, 161.950574, 322.0);
    static const Species Sm_162("Sm-162", "Sm", 38, 100, 62, 162, 8107.5745, "B-", 4343.8905, 161.941621687, 3.782);
    static const Species Eu_162("Eu-162", "Eu", 36, 99, 63, 162, 8129.5593, "B-", 5557.7761, 161.936958329, 1.41);
    static const Species Gd_162("Gd-162", "Gd", 34, 98, 64, 162, 8159.0373, "B-", 1598.8278, 161.930991812, 4.254);
    static const Species Tb_162("Tb-162", "Tb", 32, 97, 65, 162, 8164.0773, "B-", 2301.6217, 161.9292754, 2.2);
    static const Species Dy_162("Dy-162", "Dy", 30, 96, 66, 162, 8173.4555, "B-", -2140.6068, 161.92680450699999, 0.746);
    static const Species Ho_162("Ho-162", "Ho", 28, 95, 67, 162, 8155.4126, "B-", 293.6478, 161.929102543, 3.33);
    static const Species Er_162("Er-162", "Er", 26, 94, 68, 162, 8152.3959, "B-", -4856.7282, 161.928787299, 0.811);
    static const Species Tm_162("Tm-162", "Tm", 24, 93, 69, 162, 8117.5868, "B-", -1656.319, 161.934001211, 27.974);
    static const Species Yb_162("Yb-162", "Yb", 22, 92, 70, 162, 8102.5333, "B-", -6989.4042, 161.935779342, 16.213);
    static const Species Lu_162("Lu-162", "Lu", 20, 91, 71, 162, 8054.5596, "B-", -3663.3333, 161.943282776, 80.554);
    static const Species Hf_162("Hf-162", "Hf", 18, 90, 72, 162, 8027.1171, "B-", -9387.0211, 161.947215526, 9.61);
    static const Species Ta_162("Ta-162", "Ta", 16, 89, 73, 162, 7964.3432, "B-", -5782.188, 161.957292907, 67.979);
    static const Species W_162("W-162", "W", 14, 88, 74, 162, 7923.8214, "B-", -11546.0, 161.963500341, 18.955);
    static const Species Re_162("Re-162", "Re", 12, 87, 75, 162, 7848.0, "B-", -7953.0, 161.975896, 215.0);
    static const Species Os_162("Os-162", "Os", 10, 86, 76, 162, 7794.0, "B-", std::numeric_limits<double>::quiet_NaN(), 161.984434, 322.0);
    static const Species Nd_163("Nd-163", "Nd", 43, 103, 60, 163, 7992.0, "B-", 8880.0, 162.963414, 537.0);
    static const Species Pm_163("Pm-163", "Pm", 41, 102, 61, 163, 8042.0, "B-", 7640.0, 162.953881, 429.0);
    static const Species Sm_163("Sm-163", "Sm", 39, 101, 62, 163, 8084.1653, "B-", 5974.2073, 162.945679085, 7.9);
    static const Species Eu_163("Eu-163", "Eu", 37, 100, 63, 163, 8116.0172, "B-", 4814.772, 162.93926551, 0.97);
    static const Species Gd_163("Gd-163", "Gd", 35, 99, 64, 163, 8140.756, "B-", 3207.1628, 162.93409664, 0.855);
    static const Species Tb_163("Tb-163", "Tb", 33, 98, 65, 163, 8155.6322, "B-", 1785.1041, 162.930653609, 4.358);
    static const Species Dy_163("Dy-163", "Dy", 31, 97, 66, 163, 8161.784, "B-", -2.8309, 162.928737221, 0.744);
    static const Species Ho_163("Ho-163", "Ho", 29, 96, 67, 163, 8156.967, "B-", -1210.6141, 162.92874026, 0.744);
    static const Species Er_163("Er-163", "Er", 27, 95, 68, 163, 8144.7403, "B-", -2439.0, 162.930039908, 4.967);
    static const Species Tm_163("Tm-163", "Tm", 25, 94, 69, 163, 8124.9774, "B-", -3434.5345, 162.932658282, 5.92);
    static const Species Yb_163("Yb-163", "Yb", 23, 93, 70, 163, 8099.1069, "B-", -4502.4636, 162.936345406, 16.215);
    static const Species Lu_163("Lu-163", "Lu", 21, 92, 71, 163, 8066.6848, "B-", -5522.0944, 162.941179, 30.0);
    static const Species Hf_163("Hf-163", "Hf", 19, 91, 72, 163, 8028.0072, "B-", -6734.6864, 162.947107211, 27.582);
    static const Species Ta_163("Ta-163", "Ta", 17, 90, 73, 163, 7981.8905, "B-", -7626.1952, 162.954337194, 40.86);
    static const Species W_163("W-163", "W", 15, 89, 74, 163, 7930.3043, "B-", -8906.1859, 162.962524251, 62.722);
    static const Species Re_163("Re-163", "Re", 13, 88, 75, 163, 7870.8655, "B-", -9666.0, 162.972085434, 19.897);
    static const Species Os_163("Os-163", "Os", 11, 87, 76, 163, 7807.0, "B-", -11026.0, 162.982462, 322.0);
    static const Species Ir_163("Ir-163", "Ir", 9, 86, 77, 163, 7734.0, "B-", std::numeric_limits<double>::quiet_NaN(), 162.994299, 429.0);
    static const Species Pm_164("Pm-164", "Pm", 42, 103, 61, 164, 8014.0, "B-", 9565.0, 163.958819, 429.0);
    static const Species Sm_164("Sm-164", "Sm", 40, 102, 62, 164, 8067.7803, "B-", 5306.8315, 163.948550061, 4.4);
    static const Species Eu_164("Eu-164", "Eu", 38, 101, 63, 164, 8095.3686, "B-", 6461.5416, 163.942852943, 2.219);
    static const Species Gd_164("Gd-164", "Gd", 36, 100, 64, 164, 8129.9978, "B-", 2411.2959, 163.935916193, 1.073);
    static const Species Tb_164("Tb-164", "Tb", 34, 99, 65, 164, 8139.9304, "B-", 3862.6653, 163.933327561, 2.0);
    static const Species Dy_164("Dy-164", "Dy", 32, 98, 66, 164, 8158.7129, "B-", -987.1315, 163.929180819, 0.746);
    static const Species Ho_164("Ho-164", "Ho", 30, 97, 67, 164, 8147.9234, "B-", 962.0559, 163.930240548, 1.492);
    static const Species Er_164("Er-164", "Er", 28, 96, 68, 164, 8149.0191, "B-", -4033.6302, 163.929207739, 0.755);
    static const Species Tm_164("Tm-164", "Tm", 26, 95, 69, 164, 8119.6534, "B-", -896.7722, 163.933538019, 26.845);
    static const Species Yb_164("Yb-164", "Yb", 24, 94, 70, 164, 8109.4149, "B-", -6369.7952, 163.934500743, 16.217);
    static const Species Lu_164("Lu-164", "Lu", 22, 93, 71, 164, 8065.8043, "B-", -2824.0194, 163.941339, 30.0);
    static const Species Hf_164("Hf-164", "Hf", 20, 92, 72, 164, 8043.8142, "B-", -8535.5511, 163.944370709, 16.975);
    static const Species Ta_164("Ta-164", "Ta", 18, 91, 73, 164, 7986.9978, "B-", -5047.25, 163.953534, 30.0);
    static const Species W_164("W-164", "W", 16, 90, 74, 164, 7951.4515, "B-", -10763.1138, 163.958952445, 10.384);
    static const Species Re_164("Re-164", "Re", 14, 89, 75, 164, 7881.0523, "B-", -7047.718, 163.970507122, 58.566);
    static const Species Os_164("Os-164", "Os", 12, 88, 76, 164, 7833.308, "B-", -12941.0, 163.978073158, 160.927);
    static const Species Ir_164("Ir-164", "Ir", 10, 87, 77, 164, 7750.0, "B-", std::numeric_limits<double>::quiet_NaN(), 163.991966, 339.0);
    static const Species Pm_165("Pm-165", "Pm", 43, 104, 61, 165, 7992.0, "B-", 8840.0, 164.96278, 537.0);
    static const Species Sm_165("Sm-165", "Sm", 41, 103, 62, 165, 8041.0, "B-", 7219.0, 164.95329, 429.0);
    static const Species Eu_165("Eu-165", "Eu", 39, 102, 63, 165, 8080.0529, "B-", 5796.6788, 164.94554007, 5.596);
    static const Species Gd_165("Gd-165", "Gd", 37, 101, 64, 165, 8110.4428, "B-", 4063.0674, 164.93931708, 1.4);
    static const Species Tb_165("Tb-165", "Tb", 35, 100, 65, 165, 8130.3259, "B-", 3023.4392, 164.934955198, 1.654);
    static const Species Dy_165("Dy-165", "Dy", 33, 99, 66, 165, 8143.9083, "B-", 1285.7287, 164.931709402, 0.748);
    static const Species Ho_165("Ho-165", "Ho", 31, 98, 67, 165, 8146.9591, "B-", -376.6648, 164.930329116, 0.844);
    static const Species Er_165("Er-165", "Er", 29, 97, 68, 165, 8139.9348, "B-", -1591.3282, 164.930733482, 0.985);
    static const Species Tm_165("Tm-165", "Tm", 27, 96, 69, 165, 8125.5489, "B-", -2634.6364, 164.932441843, 1.779);
    static const Species Yb_165("Yb-165", "Yb", 25, 95, 70, 165, 8104.8399, "B-", -3853.1403, 164.935270241, 28.49);
    static const Species Lu_165("Lu-165", "Lu", 23, 94, 71, 165, 8076.746, "B-", -4806.735, 164.939406758, 28.49);
    static const Species Hf_165("Hf-165", "Hf", 21, 93, 72, 165, 8042.8728, "B-", -5787.641, 164.944567, 30.0);
    static const Species Ta_165("Ta-165", "Ta", 19, 92, 73, 165, 8003.0547, "B-", -6986.5555, 164.950780287, 14.571);
    static const Species W_165("W-165", "W", 17, 91, 74, 165, 7955.9704, "B-", -8201.9627, 164.958280663, 27.649);
    static const Species Re_165("Re-165", "Re", 15, 90, 75, 165, 7901.5201, "B-", -8913.0, 164.967085831, 25.328);
    static const Species Os_165("Os-165", "Os", 13, 89, 76, 165, 7843.0, "B-", -10151.0, 164.976654, 215.0);
    static const Species Ir_165("Ir-165", "Ir", 11, 88, 77, 165, 7776.0, "B-", -11277.0, 164.987552, 170.0);
    static const Species Pt_165("Pt-165", "Pt", 9, 87, 78, 165, 7703.0, "B-", std::numeric_limits<double>::quiet_NaN(), 164.999658, 429.0);
    static const Species Sm_166("Sm-166", "Sm", 42, 104, 62, 166, 8023.0, "B-", 6299.0, 165.956575, 429.0);
    static const Species Eu_166("Eu-166", "Eu", 40, 103, 63, 166, 8056.0, "B-", 7622.0, 165.949813, 107.0);
    static const Species Gd_166("Gd-166", "Gd", 38, 102, 64, 166, 8097.226, "B-", 3437.8515, 165.941630413, 1.7);
    static const Species Tb_166("Tb-166", "Tb", 36, 101, 65, 166, 8113.223, "B-", 4775.693, 165.937939727, 1.57);
    static const Species Dy_166("Dy-166", "Dy", 34, 100, 66, 166, 8137.2793, "B-", 485.8684, 165.93281281, 0.862);
    static const Species Ho_166("Ho-166", "Ho", 32, 99, 67, 166, 8135.4933, "B-", 1853.8057, 165.932291209, 0.844);
    static const Species Er_166("Er-166", "Er", 30, 98, 68, 166, 8141.9479, "B-", -3037.6667, 165.930301067, 0.358);
    static const Species Tm_166("Tm-166", "Tm", 28, 97, 69, 166, 8118.9357, "B-", -292.7714, 165.933562136, 12.401);
    static const Species Yb_166("Yb-166", "Yb", 26, 96, 70, 166, 8112.4591, "B-", -5572.7197, 165.933876439, 7.515);
    static const Species Lu_166("Lu-166", "Lu", 24, 95, 71, 166, 8074.1756, "B-", -2161.9978, 165.939859, 32.0);
    static const Species Hf_166("Hf-166", "Hf", 22, 94, 72, 166, 8056.4386, "B-", -7761.2089, 165.94218, 30.0);
    static const Species Ta_166("Ta-166", "Ta", 20, 93, 73, 166, 8004.9714, "B-", -4210.3092, 165.950512, 30.0);
    static const Species W_166("W-166", "W", 18, 92, 74, 166, 7974.8951, "B-", -10050.1358, 165.955031952, 10.159);
    static const Species Re_166("Re-166", "Re", 16, 91, 75, 166, 7909.6392, "B-", -6405.809, 165.965821216, 94.731);
    static const Species Os_166("Os-166", "Os", 14, 90, 76, 166, 7866.3371, "B-", -12126.0, 165.972698135, 19.287);
    static const Species Ir_166("Ir-166", "Ir", 12, 89, 77, 166, 7789.0, "B-", -8523.0, 165.985716, 215.0);
    static const Species Pt_166("Pt-166", "Pt", 10, 88, 78, 166, 7733.0, "B-", std::numeric_limits<double>::quiet_NaN(), 165.994866, 322.0);
    static const Species Sm_167("Sm-167", "Sm", 43, 105, 62, 167, 7992.0, "B-", 8440.0, 166.962072, 537.0);
    static const Species Eu_167("Eu-167", "Eu", 41, 104, 63, 167, 8038.0, "B-", 7006.0, 166.953011, 429.0);
    static const Species Gd_167("Gd-167", "Gd", 39, 103, 64, 167, 8075.5428, "B-", 5107.3505, 166.945490012, 5.596);
    static const Species Tb_167("Tb-167", "Tb", 37, 102, 65, 167, 8101.441, "B-", 4028.3686, 166.940007046, 2.071);
    static const Species Dy_167("Dy-167", "Dy", 35, 101, 66, 167, 8120.8782, "B-", 2368.0076, 166.935682415, 4.3);
    static const Species Ho_167("Ho-167", "Ho", 33, 100, 67, 167, 8130.3732, "B-", 1009.7971, 166.933140254, 5.57);
    static const Species Er_167("Er-167", "Er", 31, 99, 68, 167, 8131.7352, "B-", -746.1392, 166.932056192, 0.306);
    static const Species Tm_167("Tm-167", "Tm", 29, 98, 69, 167, 8122.5826, "B-", -1953.2162, 166.932857206, 1.35);
    static const Species Yb_167("Yb-167", "Yb", 27, 97, 70, 167, 8106.202, "B-", -3063.619, 166.934954069, 4.251);
    static const Species Lu_167("Lu-167", "Lu", 25, 96, 71, 167, 8083.1722, "B-", -4058.5198, 166.938243, 40.0);
    static const Species Hf_167("Hf-167", "Hf", 23, 95, 72, 167, 8054.185, "B-", -5116.6971, 166.9426, 30.0);
    static const Species Ta_167("Ta-167", "Ta", 21, 94, 73, 167, 8018.8614, "B-", -6257.8523, 166.948093, 30.0);
    static const Species W_167("W-167", "W", 19, 93, 74, 167, 7976.7045, "B-", -7259.0, 166.95481108, 20.078);
    static const Species Re_167("Re-167", "Re", 17, 92, 75, 167, 7929.0, "B-", -8335.0, 166.962604, 43.0);
    static const Species Os_167("Os-167", "Os", 15, 91, 76, 167, 7873.9557, "B-", -9426.4121, 166.971552304, 86.841);
    static const Species Ir_167("Ir-167", "Ir", 13, 90, 77, 167, 7812.8254, "B-", -10319.0, 166.981671973, 19.694);
    static const Species Pt_167("Pt-167", "Pt", 11, 89, 78, 167, 7746.0, "B-", std::numeric_limits<double>::quiet_NaN(), 166.99275, 329.0);
    static const Species Sm_168("Sm-168", "Sm", 44, 106, 62, 168, 7971.0, "B-", 7610.0, 167.966033, 322.0);
    static const Species Eu_168("Eu-168", "Eu", 42, 105, 63, 168, 8012.0, "B-", 8899.0, 167.957863, 429.0);
    static const Species Gd_168("Gd-168", "Gd", 40, 104, 64, 168, 8060.0, "B-", 4631.0, 167.948309, 322.0);
    static const Species Tb_168("Tb-168", "Tb", 38, 103, 65, 168, 8082.798, "B-", 5777.2167, 167.943337074, 4.5);
    static const Species Dy_168("Dy-168", "Dy", 36, 102, 66, 168, 8112.5293, "B-", 1500.8333, 167.937134977, 150.303);
    static const Species Ho_168("Ho-168", "Ho", 34, 101, 67, 168, 8116.8061, "B-", 2930.0, 167.935523766, 32.207);
    static const Species Er_168("Er-168", "Er", 32, 100, 68, 168, 8129.5897, "B-", -1676.8526, 167.932378282, 0.28);
    static const Species Tm_168("Tm-168", "Tm", 30, 99, 69, 168, 8114.9516, "B-", 267.488, 167.934178457, 1.8);
    static const Species Yb_168("Yb-168", "Yb", 28, 98, 70, 168, 8111.887, "B-", -4507.0351, 167.933891297, 0.1);
    static const Species Lu_168("Lu-168", "Lu", 26, 97, 71, 168, 8080.4026, "B-", -1712.274, 167.938729798, 40.766);
    static const Species Hf_168("Hf-168", "Hf", 24, 96, 72, 168, 8065.5536, "B-", -6966.6444, 167.940568, 30.0);
    static const Species Ta_168("Ta-168", "Ta", 22, 95, 73, 168, 8019.4287, "B-", -3500.9828, 167.948047, 30.0);
    static const Species W_168("W-168", "W", 20, 94, 74, 168, 7993.9327, "B-", -9098.0411, 167.951805459, 14.233);
    static const Species Re_168("Re-168", "Re", 18, 93, 75, 168, 7935.1208, "B-", -5799.8944, 167.961572607, 33.087);
    static const Species Os_168("Os-168", "Os", 16, 92, 76, 168, 7895.9408, "B-", -11328.7644, 167.96779905, 10.631);
    static const Species Ir_168("Ir-168", "Ir", 14, 91, 77, 168, 7823.8509, "B-", -7656.1526, 167.979960978, 59.277);
    static const Species Pt_168("Pt-168", "Pt", 12, 90, 78, 168, 7773.6217, "B-", -13540.0, 167.988180196, 160.96);
    static const Species Au_168("Au-168", "Au", 10, 89, 79, 168, 7688.0, "B-", std::numeric_limits<double>::quiet_NaN(), 168.002716, 429.0);
    static const Species Eu_169("Eu-169", "Eu", 43, 106, 63, 169, 7991.0, "B-", 8230.0, 168.961717, 537.0);
    static const Species Gd_169("Gd-169", "Gd", 41, 105, 64, 169, 8035.0, "B-", 6590.0, 168.952882, 429.0);
    static const Species Tb_169("Tb-169", "Tb", 39, 104, 65, 169, 8069.0, "B-", 5116.0, 168.945807, 322.0);
    static const Species Dy_169("Dy-169", "Dy", 37, 103, 66, 169, 8094.7566, "B-", 3200.0, 168.940315231, 322.781);
    static const Species Ho_169("Ho-169", "Ho", 35, 102, 67, 169, 8109.0622, "B-", 2125.1534, 168.93687989, 21.522);
    static const Species Er_169("Er-169", "Er", 33, 101, 68, 169, 8117.0078, "B-", 353.491, 168.934598444, 0.326);
    static const Species Tm_169("Tm-169", "Tm", 31, 100, 69, 169, 8114.4702, "B-", -899.127, 168.934218956, 0.792);
    static const Species Yb_169("Yb-169", "Yb", 29, 99, 70, 169, 8104.5206, "B-", -2293.0, 168.935184208, 0.191);
    static const Species Lu_169("Lu-169", "Lu", 27, 98, 71, 169, 8086.3233, "B-", -3365.6321, 168.937645845, 3.226);
    static const Species Hf_169("Hf-169", "Hf", 25, 97, 72, 169, 8061.7791, "B-", -4426.46, 168.941259, 30.0);
    static const Species Ta_169("Ta-169", "Ta", 23, 96, 73, 169, 8030.9577, "B-", -5372.5686, 168.946011, 30.0);
    static const Species W_169("W-169", "W", 21, 95, 74, 169, 7994.5381, "B-", -6508.6195, 168.951778689, 16.571);
    static const Species Re_169("Re-169", "Re", 19, 94, 75, 169, 7951.3963, "B-", -7686.2626, 168.958765979, 12.204);
    static const Species Os_169("Os-169", "Os", 17, 93, 76, 169, 7901.2862, "B-", -8629.5684, 168.967017521, 27.847);
    static const Species Ir_169("Ir-169", "Ir", 15, 92, 77, 169, 7845.5944, "B-", -9629.0, 168.976281743, 25.02);
    static const Species Pt_169("Pt-169", "Pt", 13, 91, 78, 169, 7784.0, "B-", -10676.0, 168.986619, 215.0);
    static const Species Au_169("Au-169", "Au", 11, 90, 79, 169, 7716.0, "B-", std::numeric_limits<double>::quiet_NaN(), 168.99808, 320.0);
    static const Species Eu_170("Eu-170", "Eu", 44, 107, 63, 170, 7963.0, "B-", 9989.0, 169.96687, 537.0);
    static const Species Gd_170("Gd-170", "Gd", 42, 106, 64, 170, 8017.0, "B-", 5860.0, 169.956146, 537.0);
    static const Species Tb_170("Tb-170", "Tb", 40, 105, 65, 170, 8047.0, "B-", 7000.0, 169.949855, 322.0);
    static const Species Dy_170("Dy-170", "Dy", 38, 104, 66, 170, 8084.0, "B-", 2528.0, 169.94234, 215.0);
    static const Species Ho_170("Ho-170", "Ho", 36, 103, 67, 170, 8093.7902, "B-", 3870.0, 169.939626548, 53.697);
    static const Species Er_170("Er-170", "Er", 34, 102, 68, 170, 8111.9529, "B-", -312.1997, 169.935471933, 1.488);
    static const Species Tm_170("Tm-170", "Tm", 32, 101, 69, 170, 8105.5144, "B-", 968.6148, 169.935807093, 0.785);
    static const Species Yb_170("Yb-170", "Yb", 30, 100, 70, 170, 8106.6101, "B-", -3457.695, 169.934767242, 0.011);
    static const Species Lu_170("Lu-170", "Lu", 28, 99, 71, 170, 8081.6686, "B-", -1052.3734, 169.93847923, 18.081);
    static const Species Hf_170("Hf-170", "Hf", 26, 98, 72, 170, 8070.8762, "B-", -6116.1903, 169.939609, 30.0);
    static const Species Ta_170("Ta-170", "Ta", 24, 97, 73, 170, 8030.2965, "B-", -2846.8656, 169.946175, 30.0);
    static const Species W_170("W-170", "W", 22, 96, 74, 170, 8008.9482, "B-", -8386.8083, 169.949231235, 14.165);
    static const Species Re_170("Re-170", "Re", 20, 95, 75, 170, 7955.012, "B-", -4978.3046, 169.958234844, 12.267);
    static const Species Os_170("Os-170", "Os", 18, 94, 76, 170, 7921.1258, "B-", -10743.0, 169.963579273, 10.476);
    static const Species Ir_170("Ir-170", "Ir", 16, 93, 77, 170, 7853.0, "B-", -6883.0, 169.975113, 109.0);
    static const Species Pt_170("Pt-170", "Pt", 14, 92, 78, 170, 7808.2365, "B-", -12596.0, 169.982502087, 19.588);
    static const Species Au_170("Au-170", "Au", 12, 91, 79, 170, 7730.0, "B-", -9119.0, 169.996024, 216.0);
    static const Species Hg_170("Hg-170", "Hg", 10, 90, 80, 170, 7671.0, "B-", std::numeric_limits<double>::quiet_NaN(), 170.005814, 324.0);
    static const Species Gd_171("Gd-171", "Gd", 43, 107, 64, 171, 7990.0, "B-", 7560.0, 170.961127, 537.0);
    static const Species Tb_171("Tb-171", "Tb", 41, 106, 65, 171, 8030.0, "B-", 6240.0, 170.953011, 429.0);
    static const Species Dy_171("Dy-171", "Dy", 39, 105, 66, 171, 8062.0, "B-", 4508.0, 170.946312, 215.0);
    static const Species Ho_171("Ho-171", "Ho", 37, 104, 67, 171, 8083.6021, "B-", 3200.0, 170.941472713, 644.128);
    static const Species Er_171("Er-171", "Er", 35, 103, 68, 171, 8097.7404, "B-", 1492.449, 170.938037372, 1.511);
    static const Species Tm_171("Tm-171", "Tm", 33, 102, 69, 171, 8101.8931, "B-", 96.5468, 170.936435162, 1.043);
    static const Species Yb_171("Yb-171", "Yb", 31, 101, 70, 171, 8097.8826, "B-", -1478.3526, 170.936331515, 0.013);
    static const Species Lu_171("Lu-171", "Lu", 29, 100, 71, 171, 8084.6621, "B-", -2397.1144, 170.937918591, 1.999);
    static const Species Hf_171("Hf-171", "Hf", 27, 99, 72, 171, 8066.0687, "B-", -3711.0725, 170.940492, 31.0);
    static const Species Ta_171("Ta-171", "Ta", 25, 98, 73, 171, 8039.7914, "B-", -4634.1832, 170.944476, 30.0);
    static const Species W_171("W-171", "W", 23, 97, 74, 171, 8008.1158, "B-", -5835.8106, 170.949451, 30.0);
    static const Species Re_171("Re-171", "Re", 21, 96, 75, 171, 7969.4131, "B-", -6953.0469, 170.955716, 30.0);
    static const Species Os_171("Os-171", "Os", 19, 95, 76, 171, 7924.1769, "B-", -7885.2079, 170.963180402, 19.589);
    static const Species Ir_171("Ir-171", "Ir", 17, 94, 77, 171, 7873.4895, "B-", -8945.4613, 170.97164552, 41.295);
    static const Species Pt_171("Pt-171", "Pt", 15, 93, 78, 171, 7816.6017, "B-", -9904.2655, 170.981248868, 86.904);
    static const Species Au_171("Au-171", "Au", 13, 92, 79, 171, 7754.1069, "B-", -10901.0, 170.991881533, 22.236);
    static const Species Hg_171("Hg-171", "Hg", 11, 91, 80, 171, 7686.0, "B-", std::numeric_limits<double>::quiet_NaN(), 171.003585, 329.0);
    static const Species Gd_172("Gd-172", "Gd", 44, 108, 64, 172, 7972.0, "B-", 6720.0, 171.964605, 322.0);
    static const Species Tb_172("Tb-172", "Tb", 42, 107, 65, 172, 8006.0, "B-", 8070.0, 171.957391, 537.0);
    static const Species Dy_172("Dy-172", "Dy", 40, 106, 66, 172, 8049.0, "B-", 3724.0, 171.948728, 322.0);
    static const Species Ho_172("Ho-172", "Ho", 38, 105, 67, 172, 8066.0, "B-", 4999.0, 171.94473, 210.0);
    static const Species Er_172("Er-172", "Er", 36, 104, 68, 172, 8090.4052, "B-", 890.9756, 171.939363461, 4.253);
    static const Species Tm_172("Tm-172", "Tm", 34, 103, 69, 172, 8091.0367, "B-", 1881.9024, 171.938406959, 5.884);
    static const Species Yb_172("Yb-172", "Yb", 32, 102, 70, 172, 8097.4295, "B-", -2519.3805, 171.936386654, 0.014);
    static const Species Lu_172("Lu-172", "Lu", 30, 101, 71, 172, 8078.2334, "B-", -333.8443, 171.93909132, 2.507);
    static const Species Hf_172("Hf-172", "Hf", 28, 100, 72, 172, 8071.7439, "B-", -5072.249, 171.939449716, 26.224);
    static const Species Ta_172("Ta-172", "Ta", 26, 99, 73, 172, 8037.7056, "B-", -2232.7914, 171.944895, 30.0);
    static const Species W_172("W-172", "W", 24, 98, 74, 172, 8020.1757, "B-", -7530.3525, 171.947292, 30.0);
    static const Species Re_172("Re-172", "Re", 22, 97, 75, 172, 7971.846, "B-", -4323.1979, 171.955376165, 38.183);
    static const Species Os_172("Os-172", "Os", 20, 96, 76, 172, 7942.1626, "B-", -9864.2673, 171.960017309, 13.704);
    static const Species Ir_172("Ir-172", "Ir", 18, 95, 77, 172, 7880.2637, "B-", -6272.7035, 171.970607035, 34.785);
    static const Species Pt_172("Pt-172", "Pt", 16, 94, 78, 172, 7839.246, "B-", -11788.6592, 171.977341059, 11.139);
    static const Species Au_172("Au-172", "Au", 14, 93, 79, 172, 7766.1587, "B-", -8256.6495, 171.989996704, 60.287);
    static const Species Hg_172("Hg-172", "Hg", 12, 92, 80, 172, 7713.6064, "B-", std::numeric_limits<double>::quiet_NaN(), 171.998860581, 161.098);
    static const Species Tb_173("Tb-173", "Tb", 43, 108, 65, 173, 7988.0, "B-", 7230.0, 172.960805, 537.0);
    static const Species Dy_173("Dy-173", "Dy", 41, 107, 66, 173, 8026.0, "B-", 5610.0, 172.953043, 429.0);
    static const Species Ho_173("Ho-173", "Ho", 39, 106, 67, 173, 8054.0, "B-", 4304.0, 172.94702, 320.0);
    static const Species Er_173("Er-173", "Er", 37, 105, 68, 173, 8074.0, "B-", 2602.0, 172.9424, 210.0);
    static const Species Tm_173("Tm-173", "Tm", 35, 104, 69, 173, 8084.4633, "B-", 1295.1669, 172.93960663, 4.723);
    static const Species Yb_173("Yb-173", "Yb", 33, 103, 70, 173, 8087.4276, "B-", -670.2201, 172.938216211, 0.012);
    static const Species Lu_173("Lu-173", "Lu", 31, 102, 71, 173, 8079.0312, "B-", -1469.2244, 172.938935722, 1.682);
    static const Species Hf_173("Hf-173", "Hf", 29, 101, 72, 173, 8066.0164, "B-", -3015.2464, 172.940513, 30.0);
    static const Species Ta_173("Ta-173", "Ta", 27, 100, 73, 173, 8044.065, "B-", -3669.1553, 172.94375, 30.0);
    static const Species W_173("W-173", "W", 25, 99, 74, 173, 8018.3337, "B-", -5173.5182, 172.947689, 30.0);
    static const Species Re_173("Re-173", "Re", 23, 98, 75, 173, 7983.9068, "B-", -6115.6196, 172.953243, 30.0);
    static const Species Os_173("Os-173", "Os", 21, 97, 76, 173, 7944.0341, "B-", -7169.7944, 172.959808387, 16.059);
    static const Species Ir_173("Ir-173", "Ir", 19, 96, 77, 173, 7898.068, "B-", -8331.6974, 172.967505477, 11.316);
    static const Species Pt_173("Pt-173", "Pt", 17, 95, 78, 173, 7845.3856, "B-", -9104.7415, 172.976449922, 68.096);
    static const Species Au_173("Au-173", "Au", 15, 94, 79, 173, 7788.2348, "B-", -10171.0, 172.986224263, 24.458);
    static const Species Hg_173("Hg-173", "Hg", 13, 93, 80, 173, 7725.0, "B-", std::numeric_limits<double>::quiet_NaN(), 172.997143, 215.0);
    static const Species Tb_174("Tb-174", "Tb", 44, 109, 65, 174, 7963.0, "B-", 9160.0, 173.965679, 537.0);
    static const Species Dy_174("Dy-174", "Dy", 42, 108, 66, 174, 8011.0, "B-", 4739.0, 173.955845, 537.0);
    static const Species Ho_174("Ho-174", "Ho", 40, 107, 67, 174, 8034.0, "B-", 6080.0, 173.950757, 322.0);
    static const Species Er_174("Er-174", "Er", 38, 106, 68, 174, 8064.0, "B-", 1915.0, 173.94423, 320.0);
    static const Species Tm_174("Tm-174", "Tm", 36, 105, 69, 174, 8070.6432, "B-", 3080.0, 173.942174061, 48.01);
    static const Species Yb_174("Yb-174", "Yb", 34, 104, 70, 174, 8083.8481, "B-", -1374.2287, 173.938867545, 0.011);
    static const Species Lu_174("Lu-174", "Lu", 32, 103, 71, 174, 8071.454, "B-", 274.2911, 173.94034284, 1.682);
    static const Species Hf_174("Hf-174", "Hf", 30, 102, 72, 174, 8068.5341, "B-", -4103.8117, 173.940048377, 2.425);
    static const Species Ta_174("Ta-174", "Ta", 28, 101, 73, 174, 8040.4528, "B-", -1513.6779, 173.944454, 30.0);
    static const Species W_174("W-174", "W", 26, 100, 74, 174, 8027.2572, "B-", -6553.9925, 173.946079, 30.0);
    static const Species Re_174("Re-174", "Re", 24, 99, 75, 174, 7985.0944, "B-", -3677.7177, 173.953115, 30.0);
    static const Species Os_174("Os-174", "Os", 22, 98, 76, 174, 7959.4618, "B-", -9209.4466, 173.957063192, 11.008);
    static const Species Ir_174("Ir-174", "Ir", 20, 97, 77, 174, 7902.0377, "B-", -5468.3293, 173.966949939, 12.046);
    static const Species Pt_174("Pt-174", "Pt", 18, 96, 78, 174, 7866.1143, "B-", -11259.0, 173.972820431, 11.098);
    static const Species Au_174("Au-174", "Au", 16, 95, 79, 174, 7797.0, "B-", -7417.0, 173.984908, 109.0);
    static const Species Hg_174("Hg-174", "Hg", 14, 94, 80, 174, 7749.7851, "B-", std::numeric_limits<double>::quiet_NaN(), 173.992870575, 20.623);
    static const Species Dy_175("Dy-175", "Dy", 43, 109, 66, 175, 7986.0, "B-", 6570.0, 174.960569, 537.0);
    static const Species Ho_175("Ho-175", "Ho", 41, 108, 67, 175, 8019.0, "B-", 5352.0, 174.953516, 429.0);
    static const Species Er_175("Er-175", "Er", 39, 107, 68, 175, 8045.0, "B-", 3659.0, 174.94777, 430.0);
    static const Species Tm_175("Tm-175", "Tm", 37, 106, 69, 175, 8061.7673, "B-", 2385.0, 174.94384231, 53.677);
    static const Species Yb_175("Yb-175", "Yb", 35, 105, 70, 175, 8070.9253, "B-", 470.1219, 174.941281907, 0.076);
    static const Species Lu_175("Lu-175", "Lu", 33, 104, 71, 175, 8069.1412, "B-", -683.9154, 174.940777211, 1.295);
    static const Species Hf_175("Hf-175", "Hf", 31, 103, 72, 175, 8060.7625, "B-", -2073.1103, 174.941511424, 2.45);
    static const Species Ta_175("Ta-175", "Ta", 29, 102, 73, 175, 8044.4456, "B-", -2775.8524, 174.943737, 30.0);
    static const Species W_175("W-175", "W", 27, 101, 74, 175, 8024.113, "B-", -4344.4885, 174.946717, 30.0);
    static const Species Re_175("Re-175", "Re", 25, 100, 75, 175, 7994.8168, "B-", -5182.9513, 174.951381, 30.0);
    static const Species Os_175("Os-175", "Os", 23, 99, 76, 175, 7960.7294, "B-", -6710.8496, 174.956945126, 12.64);
    static const Species Ir_175("Ir-175", "Ir", 21, 98, 77, 175, 7917.9112, "B-", -7685.8265, 174.964149519, 13.295);
    static const Species Pt_175("Pt-175", "Pt", 19, 97, 78, 175, 7869.5216, "B-", -8304.998, 174.972400593, 19.982);
    static const Species Au_175("Au-175", "Au", 17, 96, 79, 175, 7817.5939, "B-", -9434.2436, 174.981316375, 41.399);
    static const Species Hg_175("Hg-175", "Hg", 15, 95, 80, 175, 7759.2134, "B-", std::numeric_limits<double>::quiet_NaN(), 174.991444451, 87.047);
    static const Species Dy_176("Dy-176", "Dy", 44, 110, 66, 176, 7969.0, "B-", 5780.0, 175.963918, 537.0);
    static const Species Ho_176("Ho-176", "Ho", 42, 109, 67, 176, 7997.0, "B-", 7241.0, 175.957713, 537.0);
    static const Species Er_176("Er-176", "Er", 40, 108, 68, 176, 8034.0, "B-", 2741.0, 175.94994, 430.0);
    static const Species Tm_176("Tm-176", "Tm", 38, 107, 69, 176, 8045.1214, "B-", 4120.0, 175.946997707, 107.354);
    static const Species Yb_176("Yb-176", "Yb", 36, 106, 70, 176, 8064.0853, "B-", -108.9895, 175.942574706, 0.015);
    static const Species Lu_176("Lu-176", "Lu", 34, 105, 71, 176, 8059.0209, "B-", 1194.0947, 175.942691711, 1.301);
    static const Species Hf_176("Hf-176", "Hf", 32, 104, 72, 176, 8061.3604, "B-", -3211.0484, 175.941409797, 1.591);
    static const Species Ta_176("Ta-176", "Ta", 30, 103, 73, 176, 8038.6706, "B-", -723.7709, 175.944857, 33.0);
    static const Species W_176("W-176", "W", 28, 102, 74, 176, 8030.1131, "B-", -5578.7182, 175.945634, 30.0);
    static const Species Re_176("Re-176", "Re", 26, 101, 75, 176, 7993.9707, "B-", -2931.7058, 175.951623, 30.0);
    static const Species Os_176("Os-176", "Os", 24, 100, 76, 176, 7972.8681, "B-", -8249.2618, 175.954770315, 11.754);
    static const Species Ir_176("Ir-176", "Ir", 22, 99, 77, 176, 7921.5522, "B-", -4948.0041, 175.963626261, 8.679);
    static const Species Pt_176("Pt-176", "Pt", 20, 98, 78, 176, 7888.9934, "B-", -10412.9516, 175.968938162, 13.647);
    static const Species Au_176("Au-176", "Au", 18, 97, 79, 176, 7825.3837, "B-", -6736.3283, 175.980116925, 35.625);
    static const Species Hg_176("Hg-176", "Hg", 16, 96, 80, 176, 7782.664, "B-", -12369.3668, 175.98734867, 11.936);
    static const Species Tl_176("Tl-176", "Tl", 14, 95, 81, 176, 7707.9383, "B-", std::numeric_limits<double>::quiet_NaN(), 176.000627731, 89.166);
    static const Species Ho_177("Ho-177", "Ho", 43, 110, 67, 177, 7980.0, "B-", 6578.0, 176.961052, 537.0);
    static const Species Er_177("Er-177", "Er", 41, 109, 68, 177, 8013.0, "B-", 4711.0, 176.95399, 540.0);
    static const Species Tm_177("Tm-177", "Tm", 39, 108, 69, 177, 8035.0, "B-", 3417.0, 176.948932, 215.0);
    static const Species Yb_177("Yb-177", "Yb", 37, 107, 70, 177, 8049.9741, "B-", 1397.4983, 176.945263846, 0.236);
    static const Species Lu_177("Lu-177", "Lu", 35, 106, 71, 177, 8053.4495, "B-", 496.8425, 176.94376357, 1.31);
    static const Species Hf_177("Hf-177", "Hf", 33, 105, 72, 177, 8051.8365, "B-", -1166.0, 176.943230187, 1.514);
    static const Species Ta_177("Ta-177", "Ta", 31, 104, 73, 177, 8040.8289, "B-", -2013.0144, 176.94448194, 3.558);
    static const Species W_177("W-177", "W", 29, 103, 74, 177, 8025.0359, "B-", -3432.5558, 176.946643, 30.0);
    static const Species Re_177("Re-177", "Re", 27, 102, 75, 177, 8001.2229, "B-", -4312.7271, 176.950328, 30.0);
    static const Species Os_177("Os-177", "Os", 25, 101, 76, 177, 7972.4371, "B-", -5909.0234, 176.954957902, 15.687);
    static const Species Ir_177("Ir-177", "Ir", 23, 100, 77, 177, 7934.6328, "B-", -6676.988, 176.9613015, 21.213);
    static const Species Pt_177("Pt-177", "Pt", 21, 99, 78, 177, 7892.4896, "B-", -7824.7003, 176.968469541, 16.09);
    static const Species Au_177("Au-177", "Au", 19, 98, 79, 177, 7843.8623, "B-", -8769.9135, 176.976869701, 10.7);
    static const Species Hg_177("Hg-177", "Hg", 17, 97, 80, 177, 7789.8947, "B-", -9435.7198, 176.98628459, 90.952);
    static const Species Tl_177("Tl-177", "Tl", 15, 96, 81, 177, 7732.1655, "B-", std::numeric_limits<double>::quiet_NaN(), 176.996414252, 23.218);
    static const Species Ho_178("Ho-178", "Ho", 44, 111, 67, 178, 7957.0, "B-", 8130.0, 177.965507, 537.0);
    static const Species Er_178("Er-178", "Er", 42, 110, 68, 178, 7999.0, "B-", 3980.0, 177.956779, 640.0);
    static const Species Tm_178("Tm-178", "Tm", 40, 109, 69, 178, 8017.0, "B-", 5437.0, 177.952506, 322.0);
    static const Species Yb_178("Yb-178", "Yb", 38, 108, 70, 178, 8042.7386, "B-", 660.7415, 177.9466694, 7.072);
    static const Species Lu_178("Lu-178", "Lu", 36, 107, 71, 178, 8042.0554, "B-", 2097.4851, 177.945960065, 2.416);
    static const Species Hf_178("Hf-178", "Hf", 34, 106, 72, 178, 8049.4438, "B-", -1837.0, 177.943708322, 1.519);
    static const Species Ta_178("Ta-178", "Ta", 32, 105, 73, 178, 8035.0, "B-", -191.0, 177.94568, 56.0);
    static const Species W_178("W-178", "W", 30, 104, 74, 178, 8029.2584, "B-", -4753.6082, 177.945885791, 16.316);
    static const Species Re_178("Re-178", "Re", 28, 103, 75, 178, 7998.1576, "B-", -2109.2143, 177.950989, 30.0);
    static const Species Os_178("Os-178", "Os", 26, 102, 76, 178, 7981.9128, "B-", -7289.9295, 177.953253334, 14.634);
    static const Species Ir_178("Ir-178", "Ir", 24, 101, 77, 178, 7936.563, "B-", -4256.8282, 177.961079395, 20.204);
    static const Species Pt_178("Pt-178", "Pt", 22, 100, 78, 178, 7908.253, "B-", -9694.4567, 177.965649288, 10.878);
    static const Species Au_178("Au-178", "Au", 20, 99, 79, 178, 7849.3946, "B-", -5987.6832, 177.976056714, 11.0);
    static const Species Hg_178("Hg-178", "Hg", 18, 98, 80, 178, 7811.3607, "B-", -11702.0, 177.982484756, 11.548);
    static const Species Tl_178("Tl-178", "Tl", 16, 97, 81, 178, 7741.0, "B-", -8187.0, 177.995047, 110.0);
    static const Species Pb_178("Pb-178", "Pb", 14, 96, 82, 178, 7690.8359, "B-", std::numeric_limits<double>::quiet_NaN(), 178.003836171, 24.889);
    static const Species Er_179("Er-179", "Er", 43, 111, 68, 179, 7976.0, "B-", 5821.0, 178.961267, 537.0);
    static const Species Tm_179("Tm-179", "Tm", 41, 110, 69, 179, 8004.0, "B-", 4739.0, 178.955018, 429.0);
    static const Species Yb_179("Yb-179", "Yb", 39, 109, 70, 179, 8026.0, "B-", 2419.0, 178.94993, 215.0);
    static const Species Lu_179("Lu-179", "Lu", 37, 108, 71, 179, 8035.0744, "B-", 1404.0231, 178.947332985, 5.528);
    static const Species Hf_179("Hf-179", "Hf", 35, 107, 72, 179, 8038.5474, "B-", -105.5801, 178.945825705, 1.52);
    static const Species Ta_179("Ta-179", "Ta", 33, 106, 73, 179, 8033.5869, "B-", -1062.2093, 178.94593905, 1.574);
    static const Species W_179("W-179", "W", 31, 105, 74, 179, 8023.2821, "B-", -2710.9347, 178.947079378, 15.644);
    static const Species Re_179("Re-179", "Re", 29, 104, 75, 179, 8003.7666, "B-", -3564.1747, 178.949989686, 26.45);
    static const Species Os_179("Os-179", "Os", 27, 103, 76, 179, 7979.4844, "B-", -4938.4182, 178.953815985, 16.645);
    static const Species Ir_179("Ir-179", "Ir", 25, 102, 77, 179, 7947.5248, "B-", -5813.592, 178.959117594, 10.489);
    static const Species Pt_179("Pt-179", "Pt", 23, 101, 78, 179, 7910.6759, "B-", -7279.5556, 178.965358742, 8.563);
    static const Species Au_179("Au-179", "Au", 21, 100, 79, 179, 7865.6374, "B-", -8055.648, 178.973173666, 12.555);
    static const Species Hg_179("Hg-179", "Hg", 19, 99, 80, 179, 7816.2631, "B-", -8663.2918, 178.981821759, 30.188);
    static const Species Tl_179("Tl-179", "Tl", 17, 98, 81, 179, 7763.4942, "B-", -10321.2401, 178.991122185, 41.495);
    static const Species Pb_179("Pb-179", "Pb", 15, 97, 82, 179, 7701.463, "B-", std::numeric_limits<double>::quiet_NaN(), 179.002202492, 87.203);
    static const Species Er_180("Er-180", "Er", 44, 112, 68, 180, 7960.0, "B-", 4990.0, 179.96438, 537.0);
    static const Species Tm_180("Tm-180", "Tm", 42, 111, 69, 180, 7983.0, "B-", 6550.0, 179.959023, 429.0);
    static const Species Yb_180("Yb-180", "Yb", 40, 110, 70, 180, 8016.0, "B-", 1956.0, 179.951991, 322.0);
    static const Species Lu_180("Lu-180", "Lu", 38, 109, 71, 180, 8022.0394, "B-", 3103.0, 179.949890744, 75.926);
    static const Species Hf_180("Hf-180", "Hf", 36, 108, 72, 180, 8034.9319, "B-", -845.8453, 179.946559537, 1.525);
    static const Species Ta_180("Ta-180", "Ta", 34, 107, 73, 180, 8025.8864, "B-", 702.6122, 179.947467589, 2.219);
    static const Species W_180("W-180", "W", 32, 106, 74, 180, 8025.4434, "B-", -3798.8793, 179.946713304, 1.545);
    static const Species Re_180("Re-180", "Re", 30, 105, 75, 180, 7999.9922, "B-", -1481.1659, 179.950791568, 22.965);
    static const Species Os_180("Os-180", "Os", 28, 104, 76, 180, 7987.4171, "B-", -6378.6679, 179.952381665, 16.878);
    static const Species Ir_180("Ir-180", "Ir", 26, 103, 77, 180, 7947.6337, "B-", -3547.6546, 179.959229446, 23.302);
    static const Species Pt_180("Pt-180", "Pt", 24, 102, 78, 180, 7923.5781, "B-", -8804.229, 179.96303801, 10.79);
    static const Species Au_180("Au-180", "Au", 22, 101, 79, 180, 7870.3194, "B-", -5375.133, 179.972489738, 5.108);
    static const Species Hg_180("Hg-180", "Hg", 20, 100, 80, 180, 7836.1111, "B-", -10860.075, 179.97826018, 13.574);
    static const Species Tl_180("Tl-180", "Tl", 18, 99, 81, 180, 7771.431, "B-", -7449.3696, 179.98991895, 75.058);
    static const Species Pb_180("Pb-180", "Pb", 16, 98, 82, 180, 7725.6993, "B-", std::numeric_limits<double>::quiet_NaN(), 179.997916177, 13.306);
    static const Species Tm_181("Tm-181", "Tm", 43, 112, 69, 181, 7969.0, "B-", 5649.0, 180.961954, 537.0);
    static const Species Yb_181("Yb-181", "Yb", 41, 111, 70, 181, 7996.0, "B-", 3709.0, 180.95589, 320.0);
    static const Species Lu_181("Lu-181", "Lu", 39, 110, 71, 181, 8011.9301, "B-", 2605.5435, 180.951908, 135.0);
    static const Species Hf_181("Hf-181", "Hf", 37, 109, 72, 181, 8022.003, "B-", 1036.1061, 180.949110834, 1.527);
    static const Species Ta_181("Ta-181", "Ta", 35, 108, 73, 181, 8023.405, "B-", -205.1193, 180.947998528, 1.692);
    static const Species W_181("W-181", "W", 33, 107, 74, 181, 8017.9494, "B-", -1716.5331, 180.948218733, 1.554);
    static const Species Re_181("Re-181", "Re", 31, 106, 75, 181, 8004.1434, "B-", -2967.4438, 180.950061507, 13.471);
    static const Species Os_181("Os-181", "Os", 29, 105, 76, 181, 7983.4263, "B-", -4086.9327, 180.953247188, 27.201);
    static const Species Ir_181("Ir-181", "Ir", 27, 104, 77, 181, 7956.5242, "B-", -5081.5379, 180.957634691, 5.631);
    static const Species Pt_181("Pt-181", "Pt", 25, 103, 78, 181, 7924.1271, "B-", -6510.3575, 180.963089946, 14.695);
    static const Species Au_181("Au-181", "Au", 23, 102, 79, 181, 7883.8359, "B-", -7210.0126, 180.970079102, 21.445);
    static const Species Hg_181("Hg-181", "Hg", 21, 101, 80, 181, 7839.6792, "B-", -7862.378, 180.977819368, 16.513);
    static const Species Tl_181("Tl-181", "Tl", 19, 100, 81, 181, 7791.9183, "B-", -9688.1182, 180.986259978, 9.771);
    static const Species Pb_181("Pb-181", "Pb", 17, 99, 82, 181, 7734.0704, "B-", std::numeric_limits<double>::quiet_NaN(), 180.9966606, 91.29);
    static const Species Tm_182("Tm-182", "Tm", 44, 113, 69, 182, 7948.0, "B-", 7410.0, 181.966194, 537.0);
    static const Species Yb_182("Yb-182", "Yb", 42, 112, 70, 182, 7984.0, "B-", 2870.0, 181.958239, 429.0);
    static const Species Lu_182("Lu-182", "Lu", 40, 111, 71, 182, 7996.0, "B-", 4280.0, 181.955158, 215.0);
    static const Species Hf_182("Hf-182", "Hf", 38, 110, 72, 182, 8014.8381, "B-", 381.0486, 181.950563684, 6.619);
    static const Species Ta_182("Ta-182", "Ta", 36, 109, 73, 182, 8012.6332, "B-", 1815.4592, 181.950154612, 1.693);
    static const Species W_182("W-182", "W", 34, 108, 74, 182, 8018.3096, "B-", -2800.0, 181.948205636, 0.799);
    static const Species Re_182("Re-182", "Re", 32, 107, 75, 182, 7998.6264, "B-", -837.0348, 181.95121156, 109.483);
    static const Species Os_182("Os-182", "Os", 30, 106, 76, 182, 7989.7287, "B-", -5557.4266, 181.952110154, 23.344);
    static const Species Ir_182("Ir-182", "Ir", 28, 105, 77, 182, 7954.8948, "B-", -2883.262, 181.958076296, 22.509);
    static const Species Pt_182("Pt-182", "Pt", 26, 104, 78, 182, 7934.7541, "B-", -7864.4447, 181.961171605, 14.057);
    static const Species Au_182("Au-182", "Au", 24, 103, 79, 182, 7887.2442, "B-", -4727.0905, 181.969614433, 20.143);
    static const Species Hg_182("Hg-182", "Hg", 22, 102, 80, 182, 7856.9726, "B-", -10249.6721, 181.974689173, 10.51);
    static const Species Tl_182("Tl-182", "Tl", 20, 101, 81, 182, 7796.3571, "B-", -6502.6567, 181.985692649, 12.856);
    static const Species Pb_182("Pb-182", "Pb", 18, 100, 82, 182, 7756.3296, "B-", std::numeric_limits<double>::quiet_NaN(), 181.992673537, 12.975);
    static const Species Yb_183("Yb-183", "Yb", 43, 113, 70, 183, 7963.0, "B-", 4716.0, 182.962426, 429.0);
    static const Species Lu_183("Lu-183", "Lu", 41, 112, 71, 183, 7984.8125, "B-", 3567.4325, 182.957363, 86.0);
    static const Species Hf_183("Hf-183", "Hf", 39, 111, 72, 183, 8000.0315, "B-", 2010.0, 182.953533203, 32.251);
    static const Species Ta_183("Ta-183", "Ta", 37, 110, 73, 183, 8006.74, "B-", 1072.1161, 182.95137538, 1.707);
    static const Species W_183("W-183", "W", 35, 109, 74, 183, 8008.3234, "B-", -556.0, 182.950224416, 0.798);
    static const Species Re_183("Re-183", "Re", 33, 108, 75, 183, 8001.0101, "B-", -2145.9028, 182.950821306, 8.625);
    static const Species Os_183("Os-183", "Os", 31, 107, 76, 183, 7985.0087, "B-", -3461.6215, 182.953125028, 53.428);
    static const Species Ir_183("Ir-183", "Ir", 29, 106, 77, 183, 7961.8176, "B-", -4428.9417, 182.956841231, 26.486);
    static const Species Pt_183("Pt-183", "Pt", 27, 105, 78, 183, 7933.3406, "B-", -5581.7092, 182.961595895, 15.261);
    static const Species Au_183("Au-183", "Au", 25, 104, 79, 183, 7898.5644, "B-", -6386.8322, 182.967588106, 10.116);
    static const Species Hg_183("Hg-183", "Hg", 23, 103, 80, 183, 7859.3885, "B-", -7217.3941, 182.974444652, 7.604);
    static const Species Tl_183("Tl-183", "Tl", 21, 102, 81, 183, 7815.6741, "B-", -9007.2534, 182.982192843, 10.017);
    static const Species Pb_183("Pb-183", "Pb", 19, 101, 82, 183, 7762.179, "B-", std::numeric_limits<double>::quiet_NaN(), 182.991862527, 31.11);
    static const Species Yb_184("Yb-184", "Yb", 44, 114, 70, 184, 7951.0, "B-", 3700.0, 183.965002, 540.0);
    static const Species Lu_184("Lu-184", "Lu", 42, 113, 71, 184, 7967.0, "B-", 5199.0, 183.96103, 215.0);
    static const Species Hf_184("Hf-184", "Hf", 40, 112, 72, 184, 7990.7228, "B-", 1340.0, 183.955448507, 42.625);
    static const Species Ta_184("Ta-184", "Ta", 38, 111, 73, 184, 7993.7535, "B-", 2866.0, 183.954009958, 27.923);
    static const Species W_184("W-184", "W", 36, 110, 74, 184, 8005.0777, "B-", -1485.6333, 183.95093318, 0.792);
    static const Species Re_184("Re-184", "Re", 34, 109, 75, 184, 7992.7517, "B-", 32.746, 183.952528073, 4.59);
    static const Species Os_184("Os-184", "Os", 32, 108, 76, 184, 7988.6778, "B-", -4641.7101, 183.952492919, 0.89);
    static const Species Ir_184("Ir-184", "Ir", 30, 107, 77, 184, 7959.1992, "B-", -2278.3688, 183.957476, 30.0);
    static const Species Pt_184("Pt-184", "Pt", 28, 106, 78, 184, 7942.5649, "B-", -7013.7724, 183.959921929, 15.828);
    static const Species Au_184("Au-184", "Au", 26, 105, 79, 184, 7900.1947, "B-", -3973.9269, 183.967451523, 23.912);
    static const Species Hg_184("Hg-184", "Hg", 24, 104, 80, 184, 7874.3454, "B-", -9461.4321, 183.971717709, 10.235);
    static const Species Tl_184("Tl-184", "Tl", 22, 103, 81, 184, 7818.6727, "B-", -5831.7688, 183.981874973, 10.747);
    static const Species Pb_184("Pb-184", "Pb", 20, 102, 82, 184, 7782.7264, "B-", -12306.0, 183.988135634, 13.743);
    static const Species Bi_184("Bi-184", "Bi", 18, 101, 83, 184, 7712.0, "B-", std::numeric_limits<double>::quiet_NaN(), 184.001347, 131.0);
    static const Species Yb_185("Yb-185", "Yb", 45, 115, 70, 185, 7929.0, "B-", 5480.0, 184.969425, 537.0);
    static const Species Lu_185("Lu-185", "Lu", 43, 114, 71, 185, 7955.0, "B-", 4359.0, 184.963542, 322.0);
    static const Species Hf_185("Hf-185", "Hf", 41, 113, 72, 185, 7973.9711, "B-", 3074.5667, 184.958862, 69.0);
    static const Species Ta_185("Ta-185", "Ta", 39, 112, 73, 185, 7986.3615, "B-", 1993.5, 184.955561317, 15.202);
    static const Species W_185("W-185", "W", 37, 111, 74, 185, 7992.9083, "B-", 431.1764, 184.953421206, 0.793);
    static const Species Re_185("Re-185", "Re", 35, 110, 75, 185, 7991.0101, "B-", -1013.1393, 184.95295832, 0.879);
    static const Species Os_185("Os-185", "Os", 33, 109, 76, 185, 7981.3047, "B-", -2470.3505, 184.954045969, 0.893);
    static const Species Ir_185("Ir-185", "Ir", 31, 108, 77, 185, 7963.7226, "B-", -3647.4137, 184.956698, 30.0);
    static const Species Pt_185("Pt-185", "Pt", 29, 107, 78, 185, 7939.7779, "B-", -4829.9942, 184.960613659, 27.731);
    static const Species Au_185("Au-185", "Au", 27, 106, 79, 185, 7909.441, "B-", -5674.4989, 184.965798871, 2.8);
    static const Species Hg_185("Hg-185", "Hg", 25, 105, 80, 185, 7874.5391, "B-", -6425.9062, 184.971890696, 14.641);
    static const Species Tl_185("Tl-185", "Tl", 23, 104, 81, 185, 7835.5756, "B-", -8216.5333, 184.978789189, 22.194);
    static const Species Pb_185("Pb-185", "Pb", 21, 103, 82, 185, 7786.933, "B-", -9305.0, 184.98761, 17.364);
    static const Species Bi_185("Bi-185", "Bi", 19, 102, 83, 185, 7732.0, "B-", std::numeric_limits<double>::quiet_NaN(), 184.9976, 87.0);
    static const Species Lu_186("Lu-186", "Lu", 44, 115, 71, 186, 7936.0, "B-", 6104.0, 185.96745, 429.0);
    static const Species Hf_186("Hf-186", "Hf", 42, 114, 72, 186, 7964.3032, "B-", 2183.3883, 185.960897, 55.0);
    static const Species Ta_186("Ta-186", "Ta", 40, 113, 73, 186, 7971.8357, "B-", 3901.0, 185.958553036, 64.425);
    static const Species W_186("W-186", "W", 38, 112, 74, 186, 7988.6026, "B-", -581.2819, 185.95436514, 1.302);
    static const Species Re_186("Re-186", "Re", 36, 111, 75, 186, 7981.2713, "B-", 1072.7114, 185.954989172, 0.88);
    static const Species Os_186("Os-186", "Os", 34, 110, 76, 186, 7982.8324, "B-", -3827.6813, 185.953837569, 0.816);
    static const Species Ir_186("Ir-186", "Ir", 32, 109, 77, 186, 7958.0473, "B-", -1307.903, 185.957946754, 17.74);
    static const Species Pt_186("Pt-186", "Pt", 30, 108, 78, 186, 7946.8094, "B-", -6149.5913, 185.959350845, 23.344);
    static const Species Au_186("Au-186", "Au", 28, 107, 79, 186, 7909.5409, "B-", -3175.7972, 185.965952703, 22.509);
    static const Species Hg_186("Hg-186", "Hg", 26, 106, 80, 186, 7888.2605, "B-", -8656.119, 185.969362061, 12.507);
    static const Species Tl_186("Tl-186", "Tl", 24, 105, 81, 186, 7837.5161, "B-", -5202.0427, 185.978654787, 22.276);
    static const Species Pb_186("Pb-186", "Pb", 22, 104, 82, 186, 7805.342, "B-", -11535.3999, 185.984239409, 11.813);
    static const Species Bi_186("Bi-186", "Bi", 20, 103, 83, 186, 7739.1175, "B-", -7247.0279, 185.996623169, 18.2);
    static const Species Po_186("Po-186", "Po", 18, 102, 84, 186, 7695.9488, "B-", std::numeric_limits<double>::quiet_NaN(), 186.004403174, 19.622);
    static const Species Lu_187("Lu-187", "Lu", 45, 116, 71, 187, 7923.0, "B-", 5230.0, 186.970188, 429.0);
    static const Species Hf_187("Hf-187", "Hf", 43, 115, 72, 187, 7947.0, "B-", 3896.0, 186.964573, 215.0);
    static const Species Ta_187("Ta-187", "Ta", 41, 114, 73, 187, 7963.2123, "B-", 3008.4937, 186.960391, 60.0);
    static const Species W_187("W-187", "W", 39, 113, 74, 187, 7975.1168, "B-", 1312.5048, 186.957161249, 1.302);
    static const Species Re_187("Re-187", "Re", 37, 112, 75, 187, 7977.9519, "B-", 2.4667, 186.955752217, 0.791);
    static const Species Os_187("Os-187", "Os", 35, 111, 76, 187, 7973.7814, "B-", -1669.6385, 186.955749569, 0.791);
    static const Species Ir_187("Ir-187", "Ir", 33, 110, 77, 187, 7960.6692, "B-", -2864.0151, 186.957542, 30.0);
    static const Species Pt_187("Pt-187", "Pt", 31, 109, 78, 187, 7941.1699, "B-", -3656.5811, 186.960616646, 25.837);
    static const Species Au_187("Au-187", "Au", 29, 108, 79, 187, 7917.4323, "B-", -4910.2713, 186.964542147, 24.153);
    static const Species Hg_187("Hg-187", "Hg", 27, 107, 80, 187, 7886.9905, "B-", -5673.917, 186.96981354, 13.81);
    static const Species Tl_187("Tl-187", "Tl", 25, 106, 81, 187, 7852.465, "B-", -7457.637, 186.97590474, 8.64);
    static const Species Pb_187("Pb-187", "Pb", 23, 105, 82, 187, 7808.401, "B-", -8603.6798, 186.983910842, 5.468);
    static const Species Bi_187("Bi-187", "Bi", 21, 104, 83, 187, 7758.2083, "B-", -9207.0832, 186.993147272, 10.74);
    static const Species Po_187("Po-187", "Po", 19, 103, 84, 187, 7704.7889, "B-", std::numeric_limits<double>::quiet_NaN(), 187.003031482, 35.03);
    static const Species Lu_188("Lu-188", "Lu", 46, 117, 71, 188, 7903.0, "B-", 7009.0, 187.974428, 429.0);
    static const Species Hf_188("Hf-188", "Hf", 44, 116, 72, 188, 7936.0, "B-", 3080.0, 187.966903, 322.0);
    static const Species Ta_188("Ta-188", "Ta", 42, 115, 73, 188, 7948.0, "B-", 4758.0, 187.963596, 215.0);
    static const Species W_188("W-188", "W", 40, 114, 74, 188, 7969.0532, "B-", 349.0, 187.958488325, 3.316);
    static const Species Re_188("Re-188", "Re", 38, 113, 75, 188, 7966.7481, "B-", 2120.4209, 187.958113658, 0.792);
    static const Species Os_188("Os-188", "Os", 36, 112, 76, 188, 7973.8656, "B-", -2792.3457, 187.955837292, 0.788);
    static const Species Ir_188("Ir-188", "Ir", 34, 111, 77, 188, 7954.8512, "B-", -523.986, 187.958834999, 10.116);
    static const Species Pt_188("Pt-188", "Pt", 32, 110, 78, 188, 7947.9027, "B-", -5449.6549, 187.959397521, 5.694);
    static const Species Au_188("Au-188", "Au", 30, 109, 79, 188, 7914.7537, "B-", -2172.9634, 187.965247966, 2.9);
    static const Species Hg_188("Hg-188", "Hg", 28, 108, 80, 188, 7899.034, "B-", -7861.9485, 187.967580738, 7.285);
    static const Species Tl_188("Tl-188", "Tl", 26, 107, 81, 188, 7853.0537, "B-", -4525.3784, 187.976020886, 32.103);
    static const Species Pb_188("Pb-188", "Pb", 24, 106, 82, 188, 7824.8211, "B-", -10616.2241, 187.980879079, 10.868);
    static const Species Bi_188("Bi-188", "Bi", 22, 105, 83, 188, 7764.1904, "B-", -6650.4226, 187.992276064, 12.001);
    static const Species Po_188("Po-188", "Po", 20, 104, 84, 188, 7724.6544, "B-", std::numeric_limits<double>::quiet_NaN(), 187.999415586, 21.438);
    static const Species Hf_189("Hf-189", "Hf", 45, 117, 72, 189, 7917.0, "B-", 4809.0, 188.970853, 322.0);
    static const Species Ta_189("Ta-189", "Ta", 43, 116, 73, 189, 7938.0, "B-", 3850.0, 188.96569, 215.0);
    static const Species W_189("W-189", "W", 41, 115, 74, 189, 7954.0, "B-", 2170.0, 188.961557, 215.0);
    static const Species Re_189("Re-189", "Re", 39, 114, 75, 189, 7961.8105, "B-", 1007.7049, 188.959227764, 8.793);
    static const Species Os_189("Os-189", "Os", 37, 113, 76, 189, 7963.0029, "B-", -537.1494, 188.958145949, 0.715);
    static const Species Ir_189("Ir-189", "Ir", 35, 112, 77, 189, 7956.0214, "B-", -1980.247, 188.958722602, 13.5);
    static const Species Pt_189("Pt-189", "Pt", 33, 111, 78, 189, 7941.4045, "B-", -2887.4471, 188.960848485, 10.832);
    static const Species Au_189("Au-189", "Au", 31, 110, 79, 189, 7921.9876, "B-", -3955.58, 188.963948286, 21.558);
    static const Species Hg_189("Hg-189", "Hg", 29, 109, 80, 189, 7896.9192, "B-", -5010.2727, 188.968194776, 33.873);
    static const Species Tl_189("Tl-189", "Tl", 27, 108, 81, 189, 7866.2704, "B-", -6772.0862, 188.973573525, 8.983);
    static const Species Pb_189("Pb-189", "Pb", 25, 107, 82, 189, 7826.2999, "B-", -7779.3555, 188.980843658, 15.096);
    static const Species Bi_189("Bi-189", "Bi", 23, 106, 83, 189, 7780.9999, "B-", -8642.6682, 188.989195139, 22.384);
    static const Species Po_189("Po-189", "Po", 21, 105, 84, 189, 7731.1321, "B-", std::numeric_limits<double>::quiet_NaN(), 188.998473425, 23.681);
    static const Species Hf_190("Hf-190", "Hf", 46, 118, 72, 190, 7905.0, "B-", 3920.0, 189.973376, 429.0);
    static const Species Ta_190("Ta-190", "Ta", 44, 117, 73, 190, 7922.0, "B-", 5649.0, 189.969168, 215.0);
    static const Species W_190("W-190", "W", 42, 116, 74, 190, 7947.5031, "B-", 1214.1824, 189.963103542, 37.993);
    static const Species Re_190("Re-190", "Re", 40, 115, 75, 190, 7949.7759, "B-", 3124.8105, 189.961800064, 5.227);
    static const Species Os_190("Os-190", "Os", 38, 114, 76, 190, 7962.1047, "B-", -1954.2108, 189.958445442, 0.697);
    static const Species Ir_190("Ir-190", "Ir", 36, 113, 77, 190, 7947.7017, "B-", 552.8893, 189.960543374, 1.47);
    static const Species Pt_190("Pt-190", "Pt", 34, 112, 78, 190, 7946.4941, "B-", -4472.9637, 189.959949823, 0.705);
    static const Species Au_190("Au-190", "Au", 32, 111, 79, 190, 7918.8345, "B-", -1462.915, 189.964751746, 3.7);
    static const Species Hg_190("Hg-190", "Hg", 30, 110, 80, 190, 7907.0174, "B-", -7004.3892, 189.96632225, 17.076);
    static const Species Tl_190("Tl-190", "Tl", 28, 109, 81, 190, 7866.0345, "B-", -3949.6288, 189.973841771, 7.784);
    static const Species Pb_190("Pb-190", "Pb", 26, 108, 82, 190, 7841.1294, "B-", -9820.7014, 189.978081872, 13.434);
    static const Species Bi_190("Bi-190", "Bi", 24, 107, 83, 190, 7785.3239, "B-", -6033.1975, 189.988624828, 22.515);
    static const Species Po_190("Po-190", "Po", 22, 106, 84, 190, 7749.4526, "B-", std::numeric_limits<double>::quiet_NaN(), 189.995101731, 14.131);
    static const Species Ta_191("Ta-191", "Ta", 45, 118, 73, 191, 7911.0, "B-", 4657.0, 190.97153, 322.0);
    static const Species W_191("W-191", "W", 43, 117, 74, 191, 7931.4359, "B-", 3174.2318, 190.966531, 45.0);
    static const Species Re_191("Re-191", "Re", 41, 116, 75, 191, 7943.9588, "B-", 2044.8311, 190.963123322, 11.019);
    static const Species Os_191("Os-191", "Os", 39, 115, 76, 191, 7950.5687, "B-", 313.5873, 190.960928105, 0.707);
    static const Species Ir_191("Ir-191", "Ir", 37, 114, 77, 191, 7948.1144, "B-", -1010.4903, 190.960591455, 1.406);
    static const Species Pt_191("Pt-191", "Pt", 35, 113, 78, 191, 7938.7279, "B-", -1900.4257, 190.961676261, 4.43);
    static const Species Au_191("Au-191", "Au", 33, 112, 79, 191, 7924.6819, "B-", -3206.0616, 190.963716452, 5.288);
    static const Species Hg_191("Hg-191", "Hg", 31, 111, 80, 191, 7903.8002, "B-", -4308.8981, 190.967158301, 23.918);
    static const Species Tl_191("Tl-191", "Tl", 29, 110, 81, 191, 7877.1445, "B-", -5991.7073, 190.971784093, 7.889);
    static const Species Pb_191("Pb-191", "Pb", 27, 109, 82, 191, 7841.6782, "B-", -7051.8922, 190.978216455, 7.099);
    static const Species Bi_191("Bi-191", "Bi", 25, 108, 83, 191, 7800.6613, "B-", -8170.6205, 190.985786972, 8.037);
    static const Species Po_191("Po-191", "Po", 23, 107, 84, 191, 7753.7871, "B-", -8932.644, 190.994558494, 7.624);
    static const Species At_191("At-191", "At", 21, 106, 85, 191, 7702.9233, "B-", std::numeric_limits<double>::quiet_NaN(), 191.004148081, 17.287);
    static const Species Ta_192("Ta-192", "Ta", 46, 119, 73, 192, 7894.0, "B-", 6520.0, 191.975201, 429.0);
    static const Species W_192("W-192", "W", 44, 118, 74, 192, 7924.0, "B-", 1969.0, 191.968202, 215.0);
    static const Species Re_192("Re-192", "Re", 42, 117, 75, 192, 7930.2389, "B-", 4293.475, 191.966088, 76.0);
    static const Species Os_192("Os-192", "Os", 40, 116, 76, 192, 7948.526, "B-", -1046.6722, 191.961478765, 2.484);
    static const Species Ir_192("Ir-192", "Ir", 38, 115, 77, 192, 7938.9999, "B-", 1452.8946, 191.962602414, 1.41);
    static const Species Pt_192("Pt-192", "Pt", 36, 114, 78, 192, 7942.4923, "B-", -3516.3415, 191.961042667, 2.758);
    static const Species Au_192("Au-192", "Au", 34, 113, 79, 192, 7920.1033, "B-", -760.7028, 191.964817615, 16.991);
    static const Species Hg_192("Hg-192", "Hg", 32, 112, 80, 192, 7912.0666, "B-", -6139.2324, 191.965634263, 16.679);
    static const Species Tl_192("Tl-192", "Tl", 30, 111, 81, 192, 7876.0167, "B-", -3320.4029, 191.972225, 34.0);
    static const Species Pb_192("Pb-192", "Pb", 28, 110, 82, 192, 7854.6482, "B-", -9017.3089, 191.975789598, 6.147);
    static const Species Bi_192("Bi-192", "Bi", 26, 109, 83, 192, 7803.6084, "B-", -5468.0534, 191.985470077, 32.326);
    static const Species Po_192("Po-192", "Po", 24, 108, 84, 192, 7771.0542, "B-", -10992.2252, 191.991340274, 11.416);
    static const Species At_192("At-192", "At", 22, 107, 85, 192, 7709.7283, "B-", std::numeric_limits<double>::quiet_NaN(), 192.003140912, 29.922);
    static const Species Ta_193("Ta-193", "Ta", 47, 120, 73, 193, 7883.0, "B-", 5380.0, 192.97766, 429.0);
    static const Species W_193("W-193", "W", 45, 119, 74, 193, 7907.0, "B-", 4042.0, 192.971884, 215.0);
    static const Species Re_193("Re-193", "Re", 43, 118, 75, 193, 7923.9378, "B-", 3162.7597, 192.967545, 42.0);
    static const Species Os_193("Os-193", "Os", 41, 117, 76, 193, 7936.2716, "B-", 1141.9038, 192.964149637, 2.49);
    static const Species Ir_193("Ir-193", "Ir", 39, 116, 77, 193, 7938.1346, "B-", -56.6276, 192.962923753, 1.425);
    static const Species Pt_193("Pt-193", "Pt", 37, 115, 78, 193, 7933.7875, "B-", -1074.8477, 192.962984546, 1.458);
    static const Species Au_193("Au-193", "Au", 35, 114, 79, 193, 7924.1648, "B-", -2342.6641, 192.964138442, 9.311);
    static const Species Hg_193("Hg-193", "Hg", 33, 113, 80, 193, 7907.973, "B-", -3584.9466, 192.966653395, 16.645);
    static const Species Tl_193("Tl-193", "Tl", 31, 112, 81, 193, 7885.3445, "B-", -5247.9637, 192.970501994, 7.2);
    static const Species Pb_193("Pb-193", "Pb", 29, 111, 82, 193, 7854.0994, "B-", -6344.6913, 192.976135914, 11.044);
    static const Species Bi_193("Bi-193", "Bi", 27, 110, 83, 193, 7817.1718, "B-", -7559.2614, 192.98294722, 8.132);
    static const Species Po_193("Po-193", "Po", 25, 109, 84, 193, 7773.951, "B-", -8257.9789, 192.991062421, 15.599);
    static const Species At_193("At-193", "At", 23, 108, 85, 193, 7727.1099, "B-", -9110.2435, 192.999927725, 23.222);
    static const Species Rn_193("Rn-193", "Rn", 21, 107, 86, 193, 7675.853, "B-", std::numeric_limits<double>::quiet_NaN(), 193.009707973, 26.958);
    static const Species Ta_194("Ta-194", "Ta", 48, 121, 73, 194, 7865.0, "B-", 7280.0, 193.98161, 537.0);
    static const Species W_194("W-194", "W", 46, 120, 74, 194, 7899.0, "B-", 2850.0, 193.973795, 322.0);
    static const Species Re_194("Re-194", "Re", 44, 119, 75, 194, 7909.0, "B-", 5175.0, 193.970735, 215.0);
    static const Species Os_194("Os-194", "Os", 42, 118, 76, 194, 7932.0232, "B-", 96.6, 193.965179407, 2.579);
    static const Species Ir_194("Ir-194", "Ir", 40, 117, 77, 194, 7928.4885, "B-", 2228.3252, 193.965075703, 1.429);
    static const Species Pt_194("Pt-194", "Pt", 38, 116, 78, 194, 7935.942, "B-", -2548.1518, 193.962683498, 0.532);
    static const Species Au_194("Au-194", "Au", 36, 115, 79, 194, 7918.7744, "B-", -27.9978, 193.965419051, 2.273);
    static const Species Hg_194("Hg-194", "Hg", 34, 114, 80, 194, 7914.5974, "B-", -5246.4542, 193.965449108, 3.1);
    static const Species Tl_194("Tl-194", "Tl", 32, 113, 81, 194, 7883.5211, "B-", -2729.6315, 193.971081408, 15.0);
    static const Species Pb_194("Pb-194", "Pb", 30, 112, 82, 194, 7865.4181, "B-", -8184.8462, 193.974011788, 18.717);
    static const Species Bi_194("Bi-194", "Bi", 28, 111, 83, 194, 7819.1955, "B-", -5018.4029, 193.982798581, 5.638);
    static const Species Po_194("Po-194", "Po", 26, 110, 84, 194, 7789.2947, "B-", -10288.1273, 193.988186058, 13.86);
    static const Species At_194("At-194", "At", 24, 109, 85, 194, 7732.2304, "B-", -6441.1134, 193.999230816, 25.23);
    static const Species Rn_194("Rn-194", "Rn", 22, 108, 86, 194, 7694.9961, "B-", std::numeric_limits<double>::quiet_NaN(), 194.006145636, 17.884);
    static const Species W_195("W-195", "W", 47, 121, 74, 195, 7881.0, "B-", 4820.0, 194.977735, 322.0);
    static const Species Re_195("Re-195", "Re", 45, 120, 75, 195, 7901.0, "B-", 3951.0, 194.97256, 322.0);
    static const Species Os_195("Os-195", "Os", 43, 119, 76, 195, 7917.7449, "B-", 2180.722, 194.968318, 60.0);
    static const Species Ir_195("Ir-195", "Ir", 41, 118, 77, 195, 7924.916, "B-", 1101.5601, 194.965976898, 1.431);
    static const Species Pt_195("Pt-195", "Pt", 39, 117, 78, 195, 7926.553, "B-", -226.8175, 194.964794325, 0.54);
    static const Species Au_195("Au-195", "Au", 37, 116, 79, 195, 7921.3778, "B-", -1553.719, 194.965037823, 1.201);
    static const Species Hg_195("Hg-195", "Hg", 35, 115, 80, 195, 7909.398, "B-", -2858.0499, 194.966705809, 24.843);
    static const Species Tl_195("Tl-195", "Tl", 33, 114, 81, 195, 7890.7293, "B-", -4417.2524, 194.969774052, 11.944);
    static const Species Pb_195("Pb-195", "Pb", 31, 113, 82, 195, 7864.0647, "B-", -5712.4729, 194.974516167, 5.461);
    static const Species Bi_195("Bi-195", "Bi", 29, 112, 83, 195, 7830.7579, "B-", -6908.9121, 194.980648759, 5.675);
    static const Species Po_195("Po-195", "Po", 27, 111, 84, 195, 7791.3155, "B-", -7646.3554, 194.988065781, 6.486);
    static const Species At_195("At-195", "At", 25, 110, 85, 195, 7748.0914, "B-", -8520.5842, 194.99627448, 10.276);
    static const Species Rn_195("Rn-195", "Rn", 23, 109, 86, 195, 7700.3841, "B-", std::numeric_limits<double>::quiet_NaN(), 195.005421703, 55.487);
    static const Species W_196("W-196", "W", 48, 122, 74, 196, 7872.0, "B-", 3620.0, 195.979882, 429.0);
    static const Species Re_196("Re-196", "Re", 46, 121, 75, 196, 7886.0, "B-", 5918.0, 195.975996, 322.0);
    static const Species Os_196("Os-196", "Os", 44, 120, 76, 196, 7912.2301, "B-", 1158.3989, 195.969643261, 43.0);
    static const Species Ir_196("Ir-196", "Ir", 42, 119, 77, 196, 7914.1487, "B-", 3209.0164, 195.968399669, 41.239);
    static const Species Pt_196("Pt-196", "Pt", 40, 118, 78, 196, 7926.5297, "B-", -1505.8204, 195.964954648, 0.547);
    static const Species Au_196("Au-196", "Au", 38, 117, 79, 196, 7914.8553, "B-", 687.2263, 195.966571213, 3.179);
    static const Species Hg_196("Hg-196", "Hg", 36, 116, 80, 196, 7914.37, "B-", -4329.3455, 195.965833445, 3.163);
    static const Species Tl_196("Tl-196", "Tl", 34, 115, 81, 196, 7888.29, "B-", -2148.3639, 195.970481189, 13.0);
    static const Species Pb_196("Pb-196", "Pb", 32, 114, 82, 196, 7873.3374, "B-", -7339.202, 195.972787552, 8.277);
    static const Species Bi_196("Bi-196", "Bi", 30, 113, 83, 196, 7831.9009, "B-", -4540.3012, 195.980666509, 26.224);
    static const Species Po_196("Po-196", "Po", 28, 112, 84, 196, 7804.7445, "B-", -9555.5564, 195.985540722, 5.778);
    static const Species At_196("At-196", "At", 26, 111, 85, 196, 7752.0001, "B-", -5888.3439, 195.995799034, 32.458);
    static const Species Rn_196("Rn-196", "Rn", 24, 110, 86, 196, 7717.966, "B-", std::numeric_limits<double>::quiet_NaN(), 196.002120431, 15.087);
    static const Species W_197("W-197", "W", 49, 123, 74, 197, 7853.0, "B-", 5480.0, 196.984036, 429.0);
    static const Species Re_197("Re-197", "Re", 47, 122, 75, 197, 7877.0, "B-", 4729.0, 196.978153, 322.0);
    static const Species Os_197("Os-197", "Os", 45, 121, 76, 197, 7897.0, "B-", 3185.0, 196.973076, 215.0);
    static const Species Ir_197("Ir-197", "Ir", 43, 120, 77, 197, 7909.0003, "B-", 2155.6519, 196.969657217, 21.588);
    static const Species Pt_197("Pt-197", "Pt", 41, 119, 78, 197, 7915.9714, "B-", 719.9769, 196.96734303, 0.575);
    static const Species Au_197("Au-197", "Au", 39, 118, 79, 197, 7915.6548, "B-", -599.5206, 196.966570103, 0.581);
    static const Species Hg_197("Hg-197", "Hg", 37, 117, 80, 197, 7908.6402, "B-", -2186.0092, 196.967213715, 3.442);
    static const Species Tl_197("Tl-197", "Tl", 35, 116, 81, 197, 7893.5725, "B-", -3608.8367, 196.969560492, 14.573);
    static const Species Pb_197("Pb-197", "Pb", 33, 115, 82, 197, 7871.2822, "B-", -5058.1894, 196.973434737, 5.157);
    static const Species Bi_197("Bi-197", "Bi", 31, 114, 83, 197, 7841.6348, "B-", -6294.1172, 196.978864927, 8.946);
    static const Species Po_197("Po-197", "Po", 29, 113, 84, 197, 7805.7136, "B-", -7037.8235, 196.985621939, 10.585);
    static const Species At_197("At-197", "At", 27, 112, 85, 197, 7766.0174, "B-", -7865.6231, 196.993177353, 8.57);
    static const Species Rn_197("Rn-197", "Rn", 25, 111, 86, 197, 7722.119, "B-", -8743.5996, 197.001621446, 17.383);
    static const Species Fr_197("Fr-197", "Fr", 23, 110, 87, 197, 7673.764, "B-", std::numeric_limits<double>::quiet_NaN(), 197.011008086, 60.584);
    static const Species Re_198("Re-198", "Re", 48, 123, 75, 198, 7861.0, "B-", 6610.0, 197.98176, 429.0);
    static const Species Os_198("Os-198", "Os", 46, 122, 76, 198, 7890.0, "B-", 2110.0, 197.974664, 215.0);
    static const Species Ir_198("Ir-198", "Ir", 44, 121, 77, 198, 7897.0, "B-", 4194.0, 197.972399, 215.0);
    static const Species Pt_198("Pt-198", "Pt", 42, 120, 78, 198, 7914.1512, "B-", -323.2251, 197.967896718, 2.254);
    static const Species Au_198("Au-198", "Au", 40, 119, 79, 198, 7908.5675, "B-", 1373.5226, 197.968243714, 0.579);
    static const Species Hg_198("Hg-198", "Hg", 38, 118, 80, 198, 7911.5532, "B-", -3425.5625, 197.966769177, 0.491);
    static const Species Tl_198("Tl-198", "Tl", 36, 117, 81, 198, 7890.3011, "B-", -1461.3103, 197.970446669, 8.1);
    static const Species Pb_198("Pb-198", "Pb", 34, 116, 82, 198, 7878.9695, "B-", -6693.5916, 197.97201545, 9.393);
    static const Species Bi_198("Bi-198", "Bi", 32, 115, 83, 198, 7841.2123, "B-", -3900.5727, 197.979201316, 29.598);
    static const Species Po_198("Po-198", "Po", 30, 114, 84, 198, 7817.5611, "B-", -8764.5316, 197.983388753, 18.705);
    static const Species At_198("At-198", "At", 28, 113, 85, 198, 7769.3446, "B-", -5478.4271, 197.992797864, 5.265);
    static const Species Rn_198("Rn-198", "Rn", 26, 112, 86, 198, 7737.7245, "B-", -10808.0177, 197.998679197, 14.406);
    static const Species Fr_198("Fr-198", "Fr", 24, 111, 87, 198, 7679.1873, "B-", std::numeric_limits<double>::quiet_NaN(), 198.010282081, 33.419);
    static const Species Re_199("Re-199", "Re", 49, 124, 75, 199, 7850.0, "B-", 5541.0, 198.984187, 429.0);
    static const Species Os_199("Os-199", "Os", 47, 123, 76, 199, 7874.0, "B-", 4128.0, 198.978239, 215.0);
    static const Species Ir_199("Ir-199", "Ir", 45, 122, 77, 199, 7891.2066, "B-", 2990.1656, 198.973807097, 44.073);
    static const Species Pt_199("Pt-199", "Pt", 43, 121, 78, 199, 7902.3011, "B-", 1705.0525, 198.970597022, 2.317);
    static const Species Au_199("Au-199", "Au", 41, 120, 79, 199, 7906.9379, "B-", 452.3142, 198.968766573, 0.581);
    static const Species Hg_199("Hg-199", "Hg", 39, 119, 80, 199, 7905.2794, "B-", -1486.6695, 198.968280994, 0.564);
    static const Species Tl_199("Tl-199", "Tl", 37, 118, 81, 199, 7893.8773, "B-", -2827.6624, 198.969877, 30.0);
    static const Species Pb_199("Pb-199", "Pb", 35, 117, 82, 199, 7875.7366, "B-", -4434.1181, 198.97291262, 7.322);
    static const Species Bi_199("Bi-199", "Bi", 33, 116, 83, 199, 7849.5232, "B-", -5558.7875, 198.977672841, 11.395);
    static const Species Po_199("Po-199", "Po", 31, 115, 84, 199, 7817.6582, "B-", -6415.4515, 198.983640445, 5.828);
    static const Species At_199("At-199", "At", 29, 114, 85, 199, 7781.4883, "B-", -7263.5308, 198.990527715, 5.78);
    static const Species Rn_199("Rn-199", "Rn", 27, 113, 86, 199, 7741.0568, "B-", -8331.2349, 198.998325436, 7.833);
    static const Species Fr_199("Fr-199", "Fr", 25, 112, 87, 199, 7695.2599, "B-", std::numeric_limits<double>::quiet_NaN(), 199.007269384, 14.734);
    static const Species Os_200("Os-200", "Os", 48, 124, 76, 200, 7867.0, "B-", 3020.0, 199.980086, 322.0);
    static const Species Ir_200("Ir-200", "Ir", 46, 123, 77, 200, 7878.0, "B-", 5030.0, 199.976844, 210.0);
    static const Species Pt_200("Pt-200", "Pt", 44, 122, 78, 200, 7899.1986, "B-", 640.9158, 199.971444609, 21.588);
    static const Species Au_200("Au-200", "Au", 42, 121, 79, 200, 7898.4915, "B-", 2263.1737, 199.970756558, 28.681);
    static const Species Hg_200("Hg-200", "Hg", 40, 120, 80, 200, 7905.8956, "B-", -2456.0403, 199.968326941, 0.568);
    static const Species Tl_200("Tl-200", "Tl", 38, 119, 81, 200, 7889.7037, "B-", -796.3695, 199.970963608, 6.182);
    static const Species Pb_200("Pb-200", "Pb", 36, 118, 82, 200, 7881.8101, "B-", -5880.2841, 199.971818546, 10.744);
    static const Species Bi_200("Bi-200", "Bi", 34, 117, 83, 200, 7848.4969, "B-", -3428.8901, 199.97813129, 24.37);
    static const Species Po_200("Po-200", "Po", 32, 116, 84, 200, 7827.4407, "B-", -7953.7897, 199.981812355, 8.136);
    static const Species At_200("At-200", "At", 30, 115, 85, 200, 7783.7601, "B-", -4987.4394, 199.990351099, 26.264);
    static const Species Rn_200("Rn-200", "Rn", 28, 114, 86, 200, 7754.9111, "B-", -10134.0327, 199.995705335, 6.217);
    static const Species Fr_200("Fr-200", "Fr", 26, 113, 87, 200, 7700.3292, "B-", std::numeric_limits<double>::quiet_NaN(), 200.006584666, 32.769);
    static const Species Os_201("Os-201", "Os", 49, 125, 76, 201, 7849.0, "B-", 5000.0, 200.984069, 322.0);
    static const Species Ir_201("Ir-201", "Ir", 47, 124, 77, 201, 7870.0, "B-", 3901.0, 200.978701, 215.0);
    static const Species Pt_201("Pt-201", "Pt", 45, 123, 78, 201, 7885.8337, "B-", 2660.0, 200.974513305, 53.788);
    static const Species Au_201("Au-201", "Au", 43, 122, 79, 201, 7895.1752, "B-", 1261.8237, 200.971657678, 3.455);
    static const Species Hg_201("Hg-201", "Hg", 41, 121, 80, 201, 7897.5607, "B-", -481.7508, 200.970303054, 0.763);
    static const Species Tl_201("Tl-201", "Tl", 39, 120, 81, 201, 7891.2717, "B-", -1909.7458, 200.970820235, 15.228);
    static const Species Pb_201("Pb-201", "Pb", 37, 119, 82, 201, 7877.8782, "B-", -3842.0269, 200.972870431, 14.758);
    static const Species Bi_201("Bi-201", "Bi", 35, 118, 83, 201, 7854.8713, "B-", -4907.8397, 200.976995017, 13.072);
    static const Species Po_201("Po-201", "Po", 33, 117, 84, 201, 7826.5619, "B-", -5731.7247, 200.982263799, 5.305);
    static const Species At_201("At-201", "At", 31, 116, 85, 201, 7794.1536, "B-", -6682.0288, 200.988417058, 8.786);
    static const Species Rn_201("Rn-201", "Rn", 29, 115, 86, 201, 7757.0174, "B-", -7695.9856, 200.995590511, 10.865);
    static const Species Fr_201("Fr-201", "Fr", 27, 114, 87, 201, 7714.8367, "B-", -8348.2439, 201.003852491, 9.747);
    static const Species Ra_201("Ra-201", "Ra", 25, 113, 88, 201, 7669.4108, "B-", std::numeric_limits<double>::quiet_NaN(), 201.012814699, 21.794);
    static const Species Os_202("Os-202", "Os", 50, 126, 76, 202, 7839.0, "B-", 4110.0, 201.986548, 429.0);
    static const Species Ir_202("Ir-202", "Ir", 48, 125, 77, 202, 7855.0, "B-", 6052.0, 201.982136, 322.0);
    static const Species Pt_202("Pt-202", "Pt", 46, 124, 78, 202, 7881.5609, "B-", 1660.854, 201.975639, 27.0);
    static const Species Au_202("Au-202", "Au", 44, 123, 79, 202, 7885.91, "B-", 2992.3278, 201.973856, 25.0);
    static const Species Hg_202("Hg-202", "Hg", 42, 122, 80, 202, 7896.8505, "B-", -1364.8906, 201.970643604, 0.757);
    static const Species Tl_202("Tl-202", "Tl", 40, 121, 81, 202, 7886.2206, "B-", -39.811, 201.972108874, 1.972);
    static const Species Pb_202("Pb-202", "Pb", 38, 120, 82, 202, 7882.1505, "B-", -5189.7539, 201.972151613, 4.075);
    static const Species Bi_202("Bi-202", "Bi", 36, 119, 83, 202, 7852.5857, "B-", -2809.285, 201.977723042, 15.032);
    static const Species Po_202("Po-202", "Po", 34, 118, 84, 202, 7834.8053, "B-", -7346.4628, 201.980738934, 9.307);
    static const Species At_202("At-202", "At", 32, 117, 85, 202, 7794.5637, "B-", -4320.5458, 201.988625686, 29.631);
    static const Species Rn_202("Rn-202", "Rn", 30, 116, 86, 202, 7769.3018, "B-", -9376.0984, 201.993263982, 18.808);
    static const Species Fr_202("Fr-202", "Fr", 28, 115, 87, 202, 7719.0125, "B-", -5973.362, 202.003329637, 6.476);
    static const Species Ra_202("Ra-202", "Ra", 26, 114, 88, 202, 7685.5684, "B-", std::numeric_limits<double>::quiet_NaN(), 202.009742305, 16.122);
    static const Species Os_203("Os-203", "Os", 51, 127, 76, 203, 7814.0, "B-", 7100.0, 202.992195, 429.0);
    static const Species Ir_203("Ir-203", "Ir", 49, 126, 77, 203, 7845.0, "B-", 5140.0, 202.984573, 429.0);
    static const Species Pt_203("Pt-203", "Pt", 47, 125, 78, 203, 7867.0, "B-", 3633.0, 202.979055, 215.0);
    static const Species Au_203("Au-203", "Au", 45, 124, 79, 203, 7880.865, "B-", 2125.7592, 202.975154492, 3.309);
    static const Species Hg_203("Hg-203", "Hg", 43, 123, 80, 203, 7887.4828, "B-", 492.1062, 202.972872396, 1.75);
    static const Species Tl_203("Tl-203", "Tl", 41, 122, 81, 203, 7886.053, "B-", -974.8265, 202.972344098, 1.257);
    static const Species Pb_203("Pb-203", "Pb", 39, 121, 82, 203, 7877.397, "B-", -3261.5896, 202.973390617, 7.036);
    static const Species Bi_203("Bi-203", "Bi", 37, 120, 83, 203, 7857.4762, "B-", -4214.0744, 202.976892077, 13.717);
    static const Species Po_203("Po-203", "Po", 35, 119, 84, 203, 7832.8632, "B-", -5148.2116, 202.981416072, 4.981);
    static const Species At_203("At-203", "At", 33, 118, 85, 203, 7803.6487, "B-", -5978.5623, 202.986942904, 11.404);
    static const Species Rn_203("Rn-203", "Rn", 31, 117, 86, 203, 7770.3437, "B-", -7060.4572, 202.993361155, 6.242);
    static const Species Fr_203("Fr-203", "Fr", 29, 116, 87, 203, 7731.7092, "B-", -7724.9182, 203.000940867, 6.689);
    static const Species Ra_203("Ra-203", "Ra", 27, 115, 88, 203, 7689.8015, "B-", std::numeric_limits<double>::quiet_NaN(), 203.009233907, 10.4);
    static const Species Ir_204("Ir-204", "Ir", 50, 127, 77, 204, 7823.0, "B-", 8050.0, 203.989726, 429.0);
    static const Species Pt_204("Pt-204", "Pt", 48, 126, 78, 204, 7859.0, "B-", 2770.0, 203.981084, 215.0);
    static const Species Au_204("Au-204", "Au", 46, 125, 79, 204, 7868.0, "B-", 4300.0, 203.97811, 215.0);
    static const Species Hg_204("Hg-204", "Hg", 44, 124, 80, 204, 7885.5455, "B-", -344.0781, 203.973494037, 0.534);
    static const Species Tl_204("Tl-204", "Tl", 42, 123, 81, 204, 7880.0238, "B-", 763.7453, 203.97386342, 1.238);
    static const Species Pb_204("Pb-204", "Pb", 40, 122, 82, 204, 7879.9326, "B-", -4463.8883, 203.973043506, 1.231);
    static const Species Bi_204("Bi-204", "Bi", 38, 121, 83, 204, 7854.2157, "B-", -2304.8814, 203.977835687, 9.854);
    static const Species Po_204("Po-204", "Po", 36, 120, 84, 204, 7839.0823, "B-", -6465.7939, 203.980310078, 10.811);
    static const Species At_204("At-204", "At", 34, 119, 85, 204, 7803.5522, "B-", -3905.136, 203.987251393, 24.335);
    static const Species Rn_204("Rn-204", "Rn", 32, 118, 86, 204, 7780.5743, "B-", -8577.4239, 203.991443729, 7.991);
    static const Species Fr_204("Fr-204", "Fr", 30, 117, 87, 204, 7734.6931, "B-", -5453.7889, 204.000651972, 26.389);
    static const Species Ra_204("Ra-204", "Ra", 28, 116, 88, 204, 7704.1238, "B-", std::numeric_limits<double>::quiet_NaN(), 204.006506855, 9.58);
    static const Species Ir_205("Ir-205", "Ir", 51, 128, 77, 205, 7805.0, "B-", 7220.0, 204.993988, 537.0);
    static const Species Pt_205("Pt-205", "Pt", 49, 127, 78, 205, 7836.0, "B-", 5750.0, 204.986237, 322.0);
    static const Species Au_205("Au-205", "Au", 47, 126, 79, 205, 7860.0, "B-", 3717.0, 204.980064, 215.0);
    static const Species Hg_205("Hg-205", "Hg", 45, 125, 80, 205, 7874.7325, "B-", 1533.0836, 204.976073151, 3.923);
    static const Species Tl_205("Tl-205", "Tl", 43, 124, 81, 205, 7878.3946, "B-", -50.6402, 204.974427318, 1.33);
    static const Species Pb_205("Pb-205", "Pb", 41, 123, 82, 205, 7874.3313, "B-", -2704.5927, 204.974481682, 1.228);
    static const Species Bi_205("Bi-205", "Bi", 39, 122, 83, 205, 7857.3218, "B-", -3544.171, 204.977385182, 5.161);
    static const Species Po_205("Po-205", "Po", 37, 121, 84, 205, 7836.2168, "B-", -4536.8793, 204.981190006, 10.798);
    static const Species At_205("At-205", "At", 35, 120, 85, 205, 7810.2694, "B-", -5274.7547, 204.986060546, 12.941);
    static const Species Rn_205("Rn-205", "Rn", 33, 119, 86, 205, 7780.7226, "B-", -6399.9479, 204.991723228, 5.453);
    static const Species Fr_205("Fr-205", "Fr", 31, 118, 87, 205, 7745.687, "B-", -7113.6693, 204.998593854, 8.399);
    static const Species Ra_205("Ra-205", "Ra", 29, 117, 88, 205, 7707.1698, "B-", -8302.8357, 205.006230692, 24.446);
    static const Species Ac_205("Ac-205", "Ac", 27, 116, 89, 205, 7662.8519, "B-", std::numeric_limits<double>::quiet_NaN(), 205.015144152, 63.682);
    static const Species Pt_206("Pt-206", "Pt", 50, 128, 78, 206, 7820.0, "B-", 4950.0, 205.99008, 322.0);
    static const Species Au_206("Au-206", "Au", 48, 127, 79, 206, 7840.0, "B-", 6755.0, 205.984766, 322.0);
    static const Species Hg_206("Hg-206", "Hg", 46, 126, 80, 206, 7869.1723, "B-", 1307.5659, 205.977513837, 21.943);
    static const Species Tl_206("Tl-206", "Tl", 44, 125, 81, 206, 7871.7219, "B-", 1532.2128, 205.976110108, 1.38);
    static const Species Pb_206("Pb-206", "Pb", 42, 124, 82, 206, 7875.362, "B-", -3757.3057, 205.97446521, 1.228);
    static const Species Bi_206("Bi-206", "Bi", 40, 123, 83, 206, 7853.3249, "B-", -1839.5323, 205.978498843, 8.193);
    static const Species Po_206("Po-206", "Po", 38, 122, 84, 206, 7840.5973, "B-", -5749.2803, 205.980473662, 4.306);
    static const Species At_206("At-206", "At", 36, 121, 85, 206, 7808.8904, "B-", -3306.4697, 205.986645768, 14.523);
    static const Species Rn_206("Rn-206", "Rn", 34, 120, 86, 206, 7789.0417, "B-", -7886.0592, 205.990195409, 9.223);
    static const Species Fr_206("Fr-206", "Fr", 32, 119, 87, 206, 7746.9621, "B-", -4812.4721, 205.998661441, 29.856);
    static const Species Ra_206("Ra-206", "Ra", 30, 118, 88, 206, 7719.8028, "B-", -9919.1406, 206.003827842, 19.332);
    static const Species Ac_206("Ac-206", "Ac", 28, 117, 89, 206, 7667.8538, "B-", std::numeric_limits<double>::quiet_NaN(), 206.014476477, 69.874);
    static const Species Pt_207("Pt-207", "Pt", 51, 129, 78, 207, 7797.0, "B-", 6501.0, 206.995556, 429.0);
    static const Species Au_207("Au-207", "Au", 49, 128, 79, 207, 7824.0, "B-", 5847.0, 206.988577, 322.0);
    static const Species Hg_207("Hg-207", "Hg", 47, 127, 80, 207, 7848.6112, "B-", 4546.9906, 206.9823, 32.0);
    static const Species Tl_207("Tl-207", "Tl", 45, 126, 81, 207, 7866.7979, "B-", 1417.5323, 206.977418605, 5.839);
    static const Species Pb_207("Pb-207", "Pb", 43, 125, 82, 207, 7869.8664, "B-", -2397.414, 206.975896821, 1.231);
    static const Species Bi_207("Bi-207", "Bi", 41, 124, 83, 207, 7854.5053, "B-", -2908.8541, 206.978470551, 2.573);
    static const Species Po_207("Po-207", "Po", 39, 123, 84, 207, 7836.6734, "B-", -3918.2186, 206.981593334, 7.148);
    static const Species At_207("At-207", "At", 37, 122, 85, 207, 7813.9653, "B-", -4592.74, 206.985799715, 13.318);
    static const Species Rn_207("Rn-207", "Rn", 35, 121, 86, 207, 7787.9987, "B-", -5785.7211, 206.990730224, 5.09);
    static const Species Fr_207("Fr-207", "Fr", 33, 120, 87, 207, 7756.2689, "B-", -6363.0084, 206.99694145, 18.847);
    static const Species Ra_207("Ra-207", "Ra", 31, 119, 88, 207, 7721.7503, "B-", -7632.2404, 207.00377242, 62.572);
    static const Species Ac_207("Ac-207", "Ac", 29, 118, 89, 207, 7681.1001, "B-", std::numeric_limits<double>::quiet_NaN(), 207.011965967, 60.384);
    static const Species Pt_208("Pt-208", "Pt", 52, 130, 78, 208, 7780.0, "B-", 5410.0, 207.999463, 429.0);
    static const Species Au_208("Au-208", "Au", 50, 129, 79, 208, 7803.0, "B-", 7355.0, 207.993655, 322.0);
    static const Species Hg_208("Hg-208", "Hg", 48, 128, 80, 208, 7834.1914, "B-", 3484.7131, 207.985759, 33.0);
    static const Species Tl_208("Tl-208", "Tl", 46, 127, 81, 208, 7847.1835, "B-", 4998.3984, 207.982018006, 1.989);
    static const Species Pb_208("Pb-208", "Pb", 44, 126, 82, 208, 7867.453, "B-", -2878.368, 207.976652005, 1.232);
    static const Species Bi_208("Bi-208", "Bi", 42, 125, 83, 208, 7849.8534, "B-", -1400.9438, 207.97974206, 2.474);
    static const Species Po_208("Po-208", "Po", 40, 124, 84, 208, 7839.3568, "B-", -4999.3061, 207.981246035, 1.795);
    static const Species At_208("At-208", "At", 38, 123, 85, 208, 7811.5604, "B-", -2814.5118, 207.986613011, 9.577);
    static const Species Rn_208("Rn-208", "Rn", 36, 122, 86, 208, 7794.2678, "B-", -6990.4614, 207.989634513, 10.91);
    static const Species Fr_208("Fr-208", "Fr", 34, 121, 87, 208, 7756.8985, "B-", -4392.8615, 207.997139082, 12.514);
    static const Species Ra_208("Ra-208", "Ra", 32, 120, 88, 208, 7732.0177, "B-", -9032.9208, 208.001855012, 9.686);
    static const Species Ac_208("Ac-208", "Ac", 30, 119, 89, 208, 7684.8289, "B-", -5927.1868, 208.011552251, 69.225);
    static const Species Th_208("Th-208", "Th", 28, 118, 90, 208, 7652.5716, "B-", std::numeric_limits<double>::quiet_NaN(), 208.017915348, 34.208);
    static const Species Au_209("Au-209", "Au", 51, 130, 79, 209, 7786.0, "B-", 6380.0, 208.997606, 429.0);
    static const Species Hg_209("Hg-209", "Hg", 49, 129, 80, 209, 7813.0, "B-", 5035.0, 208.990757, 161.0);
    static const Species Tl_209("Tl-209", "Tl", 47, 128, 81, 209, 7833.3979, "B-", 3969.7809, 208.985351713, 6.559);
    static const Species Pb_209("Pb-209", "Pb", 45, 127, 82, 209, 7848.6488, "B-", 644.0152, 208.981089978, 1.875);
    static const Species Bi_209("Bi-209", "Bi", 43, 126, 83, 209, 7847.9869, "B-", -1892.5741, 208.980398599, 1.465);
    static const Species Po_209("Po-209", "Po", 41, 125, 84, 209, 7835.1882, "B-", -3482.2417, 208.982430361, 1.909);
    static const Species At_209("At-209", "At", 39, 124, 85, 209, 7814.7835, "B-", -3942.7237, 208.986168701, 5.094);
    static const Species Rn_209("Rn-209", "Rn", 37, 123, 86, 209, 7792.1755, "B-", -5158.9051, 208.990401389, 10.692);
    static const Species Fr_209("Fr-209", "Fr", 35, 122, 87, 209, 7763.7485, "B-", -5640.3845, 208.995939701, 12.349);
    static const Species Ra_209("Ra-209", "Ra", 33, 121, 88, 209, 7733.0177, "B-", -6986.6464, 209.001994902, 6.169);
    static const Species Ac_209("Ac-209", "Ac", 31, 120, 89, 209, 7695.8455, "B-", -7550.0, 209.009495375, 59.953);
    static const Species Th_209("Th-209", "Th", 29, 119, 90, 209, 7656.0, "B-", std::numeric_limits<double>::quiet_NaN(), 209.017601, 111.0);
    static const Species Au_210("Au-210", "Au", 52, 131, 79, 210, 7764.0, "B-", 7980.0, 210.002877, 429.0);
    static const Species Hg_210("Hg-210", "Hg", 50, 130, 80, 210, 7799.0, "B-", 3947.0, 209.99431, 215.0);
    static const Species Tl_210("Tl-210", "Tl", 48, 129, 81, 210, 7813.589, "B-", 5481.4334, 209.990072942, 12.456);
    static const Species Pb_210("Pb-210", "Pb", 46, 128, 82, 210, 7835.9656, "B-", 63.4758, 209.984188381, 1.554);
    static const Species Bi_210("Bi-210", "Bi", 44, 127, 83, 210, 7832.5424, "B-", 1161.1549, 209.984120237, 1.463);
    static const Species Po_210("Po-210", "Po", 42, 126, 84, 210, 7834.3462, "B-", -3980.9605, 209.982873686, 1.23);
    static const Species At_210("At-210", "At", 40, 125, 85, 210, 7811.6638, "B-", -2367.3352, 209.987147423, 8.261);
    static const Species Rn_210("Rn-210", "Rn", 38, 124, 86, 210, 7796.6653, "B-", -6261.2558, 209.989688862, 4.892);
    static const Species Fr_210("Fr-210", "Fr", 36, 123, 87, 210, 7763.1243, "B-", -3786.3467, 209.996410596, 14.407);
    static const Species Ra_210("Ra-210", "Ra", 34, 122, 88, 210, 7741.3687, "B-", -8321.2403, 210.000475406, 9.868);
    static const Species Ac_210("Ac-210", "Ac", 32, 121, 89, 210, 7698.0182, "B-", -5295.442, 210.009408625, 66.782);
    static const Species Th_210("Th-210", "Th", 30, 120, 90, 210, 7669.0764, "B-", std::numeric_limits<double>::quiet_NaN(), 210.015093515, 20.299);
    static const Species Hg_211("Hg-211", "Hg", 51, 131, 80, 211, 7777.0, "B-", 5688.0, 210.999581, 215.0);
    static const Species Tl_211("Tl-211", "Tl", 49, 130, 81, 211, 7799.7915, "B-", 4415.0129, 210.993475, 45.0);
    static const Species Pb_211("Pb-211", "Pb", 47, 129, 82, 211, 7817.0079, "B-", 1366.1041, 210.988735288, 2.426);
    static const Species Bi_211("Bi-211", "Bi", 45, 128, 83, 211, 7819.7745, "B-", 573.3763, 210.987268715, 5.842);
    static const Species Po_211("Po-211", "Po", 43, 127, 84, 211, 7818.7842, "B-", -785.3012, 210.986653171, 1.347);
    static const Species At_211("At-211", "At", 41, 126, 85, 211, 7811.3545, "B-", -2891.8615, 210.987496226, 2.929);
    static const Species Rn_211("Rn-211", "Rn", 39, 125, 86, 211, 7793.9412, "B-", -4615.0152, 210.990600767, 7.314);
    static const Species Fr_211("Fr-211", "Fr", 37, 124, 87, 211, 7768.3613, "B-", -4972.1844, 210.995555189, 12.872);
    static const Species Ra_211("Ra-211", "Ra", 35, 123, 88, 211, 7741.0887, "B-", -6311.6152, 211.000893049, 5.331);
    static const Species Ac_211("Ac-211", "Ac", 33, 122, 89, 211, 7707.468, "B-", -6732.9111, 211.007668846, 57.706);
    static const Species Th_211("Th-211", "Th", 31, 121, 90, 211, 7671.8506, "B-", -8175.8296, 211.014896923, 92.399);
    static const Species Pa_211("Pa-211", "Pa", 29, 120, 91, 211, 7629.3948, "B-", std::numeric_limits<double>::quiet_NaN(), 211.023674036, 74.581);
    static const Species Hg_212("Hg-212", "Hg", 52, 132, 80, 212, 7762.0, "B-", 4571.0, 212.003242, 322.0);
    static const Species Tl_212("Tl-212", "Tl", 50, 131, 81, 212, 7780.0, "B-", 5998.0, 211.998335, 215.0);
    static const Species Pb_212("Pb-212", "Pb", 48, 130, 82, 212, 7804.3203, "B-", 569.0133, 211.991895891, 1.975);
    static const Species Bi_212("Bi-212", "Bi", 46, 129, 83, 212, 7803.314, "B-", 2251.4656, 211.99128503, 1.989);
    static const Species Po_212("Po-212", "Po", 44, 128, 84, 212, 7810.2438, "B-", -1741.2596, 211.988867982, 1.237);
    static const Species At_212("At-212", "At", 42, 127, 85, 212, 7798.34, "B-", 31.0705, 211.990737301, 2.559);
    static const Species Rn_212("Rn-212", "Rn", 40, 126, 86, 212, 7794.7963, "B-", -5143.221, 211.990703946, 3.338);
    static const Species Fr_212("Fr-212", "Fr", 38, 125, 87, 212, 7766.8455, "B-", -3317.2355, 211.99622542, 9.419);
    static const Species Ra_212("Ra-212", "Ra", 36, 124, 88, 212, 7747.5078, "B-", -7498.3626, 211.999786619, 11.007);
    static const Species Ac_212("Ac-212", "Ac", 34, 123, 89, 212, 7708.4479, "B-", -4811.2864, 212.007836442, 23.492);
    static const Species Th_212("Th-212", "Th", 32, 122, 90, 212, 7682.0628, "B-", -9485.6366, 212.01300157, 10.852);
    static const Species Pa_212("Pa-212", "Pa", 30, 121, 91, 212, 7633.6289, "B-", std::numeric_limits<double>::quiet_NaN(), 212.023184819, 94.047);
    static const Species Hg_213("Hg-213", "Hg", 53, 133, 80, 213, 7739.0, "B-", 6416.0, 213.008803, 322.0);
    static const Species Tl_213("Tl-213", "Tl", 51, 132, 81, 213, 7765.4311, "B-", 4987.4088, 213.001915, 29.0);
    static const Species Pb_213("Pb-213", "Pb", 49, 131, 82, 213, 7785.1732, "B-", 2028.073, 212.996560796, 7.465);
    static const Species Bi_213("Bi-213", "Bi", 47, 130, 83, 213, 7791.0217, "B-", 1421.8481, 212.99438357, 5.455);
    static const Species Po_213("Po-213", "Po", 45, 129, 84, 213, 7794.024, "B-", -73.9972, 212.992857154, 3.277);
    static const Species At_213("At-213", "At", 43, 128, 85, 213, 7790.0036, "B-", -883.5727, 212.992936593, 5.258);
    static const Species Rn_213("Rn-213", "Rn", 41, 127, 86, 213, 7782.1824, "B-", -2141.7493, 212.993885147, 3.618);
    static const Species Fr_213("Fr-213", "Fr", 39, 126, 87, 213, 7768.4543, "B-", -3899.7568, 212.99618441, 5.053);
    static const Species Ra_213("Ra-213", "Ra", 37, 125, 88, 213, 7746.4726, "B-", -5795.4721, 213.000370971, 10.54);
    static const Species Ac_213("Ac-213", "Ac", 35, 124, 89, 213, 7715.5908, "B-", -5979.0781, 213.006592665, 12.522);
    static const Species Th_213("Th-213", "Th", 33, 123, 90, 213, 7683.847, "B-", -7534.0866, 213.01301147, 9.895);
    static const Species Pa_213("Pa-213", "Pa", 31, 122, 91, 213, 7644.8027, "B-", std::numeric_limits<double>::quiet_NaN(), 213.021099644, 61.374);
    static const Species Hg_214("Hg-214", "Hg", 54, 134, 80, 214, 7724.0, "B-", 5306.0, 214.012636, 429.0);
    static const Species Tl_214("Tl-214", "Tl", 52, 133, 81, 214, 7745.0, "B-", 6648.0, 214.00694, 210.0);
    static const Species Pb_214("Pb-214", "Pb", 50, 132, 82, 214, 7772.3955, "B-", 1017.7611, 213.999803521, 2.114);
    static const Species Bi_214("Bi-214", "Bi", 48, 131, 83, 214, 7773.4955, "B-", 3269.1925, 213.998710909, 12.033);
    static const Species Po_214("Po-214", "Po", 46, 130, 84, 214, 7785.1163, "B-", -1090.8208, 213.995201287, 1.556);
    static const Species At_214("At-214", "At", 44, 129, 85, 214, 7776.3632, "B-", 940.5125, 213.996372331, 4.274);
    static const Species Rn_214("Rn-214", "Rn", 42, 128, 86, 214, 7777.1023, "B-", -3361.3369, 213.99536265, 9.862);
    static const Species Fr_214("Fr-214", "Fr", 40, 127, 87, 214, 7757.7393, "B-", -1051.0675, 213.998971193, 9.145);
    static const Species Ra_214("Ra-214", "Ra", 38, 126, 88, 214, 7749.1719, "B-", -6340.5313, 214.00009956, 5.636);
    static const Species Ac_214("Ac-214", "Ac", 36, 125, 89, 214, 7715.8874, "B-", -4261.6599, 214.0069064, 14.547);
    static const Species Th_214("Th-214", "Th", 34, 124, 90, 214, 7692.3173, "B-", -8764.963, 214.01148148, 11.445);
    static const Species Pa_214("Pa-214", "Pa", 32, 123, 91, 214, 7647.7037, "B-", std::numeric_limits<double>::quiet_NaN(), 214.020891055, 87.18);
    static const Species Hg_215("Hg-215", "Hg", 55, 135, 80, 215, 7701.0, "B-", 7079.0, 215.018368, 429.0);
    static const Species Tl_215("Tl-215", "Tl", 53, 134, 81, 215, 7730.0, "B-", 5688.0, 215.010768, 322.0);
    static const Species Pb_215("Pb-215", "Pb", 51, 133, 82, 215, 7752.7381, "B-", 2712.9729, 215.004661591, 56.56);
    static const Species Bi_215("Bi-215", "Bi", 49, 132, 83, 215, 7761.7177, "B-", 2171.0426, 215.001749095, 6.037);
    static const Species Po_215("Po-215", "Po", 47, 131, 84, 215, 7768.1768, "B-", 714.8128, 214.999418385, 2.276);
    static const Species At_215("At-215", "At", 45, 130, 85, 215, 7767.8627, "B-", -87.5935, 214.998651002, 7.116);
    static const Species Rn_215("Rn-215", "Rn", 43, 129, 86, 215, 7763.8164, "B-", -1487.1274, 214.998745037, 6.538);
    static const Species Fr_215("Fr-215", "Fr", 41, 128, 87, 215, 7753.2607, "B-", -2213.8573, 215.000341534, 7.585);
    static const Species Ra_215("Ra-215", "Ra", 39, 127, 88, 215, 7739.3249, "B-", -3498.5554, 215.002718208, 7.73);
    static const Species Ac_215("Ac-215", "Ac", 37, 126, 89, 215, 7719.4137, "B-", -4890.8833, 215.006474061, 13.318);
    static const Species Th_215("Th-215", "Th", 35, 125, 90, 215, 7693.0266, "B-", -6883.1037, 215.01172464, 6.8);
    static const Species Pa_215("Pa-215", "Pa", 33, 124, 91, 215, 7657.3733, "B-", -7084.7747, 215.019113955, 88.513);
    static const Species U_215("U-215", "U", 31, 123, 92, 215, 7620.7821, "B-", std::numeric_limits<double>::quiet_NaN(), 215.026719774, 111.794);
    static const Species Hg_216("Hg-216", "Hg", 56, 136, 80, 216, 7685.0, "B-", 6050.0, 216.022459, 429.0);
    static const Species Tl_216("Tl-216", "Tl", 54, 135, 81, 216, 7709.0, "B-", 7361.0, 216.015964, 322.0);
    static const Species Pb_216("Pb-216", "Pb", 52, 134, 82, 216, 7740.0, "B-", 1636.0, 216.008062, 215.0);
    static const Species Bi_216("Bi-216", "Bi", 50, 133, 83, 216, 7743.4996, "B-", 4091.652, 216.006305985, 12.0);
    static const Species Po_216("Po-216", "Po", 48, 132, 84, 216, 7758.8205, "B-", -474.3423, 216.001913416, 1.948);
    static const Species At_216("At-216", "At", 46, 131, 85, 216, 7753.0024, "B-", 2003.3657, 216.002422643, 3.837);
    static const Species Rn_216("Rn-216", "Rn", 44, 130, 86, 216, 7758.6553, "B-", -2717.7096, 216.000271942, 6.192);
    static const Species Fr_216("Fr-216", "Fr", 42, 129, 87, 216, 7742.4513, "B-", -320.4441, 216.003189523, 4.48);
    static const Species Ra_216("Ra-216", "Ra", 40, 128, 88, 216, 7737.3458, "B-", -4858.2701, 216.003533534, 8.592);
    static const Species Ac_216("Ac-216", "Ac", 38, 127, 89, 216, 7711.2319, "B-", -2148.8006, 216.008749101, 9.908);
    static const Species Th_216("Th-216", "Th", 36, 126, 90, 216, 7697.6617, "B-", -7525.2616, 216.011055933, 11.92);
    static const Species Pa_216("Pa-216", "Pa", 34, 125, 91, 216, 7659.2006, "B-", -5242.6308, 216.019134633, 26.459);
    static const Species U_216("U-216", "U", 32, 124, 92, 216, 7631.3072, "B-", std::numeric_limits<double>::quiet_NaN(), 216.024762829, 30.158);
    static const Species Tl_217("Tl-217", "Tl", 55, 136, 81, 217, 7693.0, "B-", 6399.0, 217.020032, 429.0);
    static const Species Pb_217("Pb-217", "Pb", 53, 135, 82, 217, 7719.0, "B-", 3530.0, 217.013162, 322.0);
    static const Species Bi_217("Bi-217", "Bi", 51, 134, 83, 217, 7731.8491, "B-", 2846.5103, 217.009372, 19.0);
    static const Species Po_217("Po-217", "Po", 49, 133, 84, 217, 7741.3614, "B-", 1488.8543, 217.006316145, 7.025);
    static const Species At_217("At-217", "At", 47, 132, 85, 217, 7744.6172, "B-", 736.032, 217.004717794, 5.368);
    static const Species Rn_217("Rn-217", "Rn", 45, 131, 86, 217, 7744.4037, "B-", -656.0967, 217.003927632, 4.506);
    static const Species Fr_217("Fr-217", "Fr", 43, 130, 87, 217, 7737.775, "B-", -1574.8729, 217.00463198, 7.011);
    static const Species Ra_217("Ra-217", "Ra", 41, 129, 88, 217, 7726.9122, "B-", -2812.7853, 217.006322676, 7.564);
    static const Species Ac_217("Ac-217", "Ac", 39, 128, 89, 217, 7710.3448, "B-", -3503.459, 217.009342325, 12.048);
    static const Species Th_217("Th-217", "Th", 37, 127, 90, 217, 7690.5945, "B-", -4848.968, 217.013103443, 11.394);
    static const Species Pa_217("Pa-217", "Pa", 35, 126, 91, 217, 7664.6438, "B-", -5916.0, 217.018309024, 13.417);
    static const Species U_217("U-217", "U", 33, 125, 92, 217, 7634.0, "B-", std::numeric_limits<double>::quiet_NaN(), 217.02466, 86.0);
    static const Species Tl_218("Tl-218", "Tl", 56, 137, 81, 218, 7672.0, "B-", 8081.0, 218.025454, 429.0);
    static const Species Pb_218("Pb-218", "Pb", 54, 136, 82, 218, 7705.0, "B-", 2414.0, 218.016779, 322.0);
    static const Species Bi_218("Bi-218", "Bi", 52, 135, 83, 218, 7712.828, "B-", 4859.3866, 218.014188, 29.0);
    static const Species Po_218("Po-218", "Po", 50, 134, 84, 218, 7731.53, "B-", 256.4334, 218.008971234, 2.112);
    static const Species At_218("At-218", "At", 48, 133, 85, 218, 7729.1175, "B-", 2882.8048, 218.008695941, 12.349);
    static const Species Rn_218("Rn-218", "Rn", 46, 132, 86, 218, 7738.7527, "B-", -1842.0267, 218.005601123, 2.486);
    static const Species Fr_218("Fr-218", "Fr", 44, 131, 87, 218, 7726.7143, "B-", 413.8838, 218.00757862, 4.546);
    static const Species Ra_218("Ra-218", "Ra", 42, 130, 88, 218, 7725.0241, "B-", -4205.2887, 218.007134297, 10.528);
    static const Species Ac_218("Ac-218", "Ac", 40, 129, 89, 218, 7702.145, "B-", -1515.9019, 218.01164886, 61.853);
    static const Species Th_218("Th-218", "Th", 38, 128, 90, 218, 7691.6026, "B-", -6282.8212, 218.013276248, 11.289);
    static const Species Pa_218("Pa-218", "Pa", 36, 127, 91, 218, 7659.1935, "B-", -3245.0869, 218.020021133, 19.158);
    static const Species U_218("U-218", "U", 34, 126, 92, 218, 7640.7191, "B-", std::numeric_limits<double>::quiet_NaN(), 218.023504877, 14.722);
    static const Species Pb_219("Pb-219", "Pb", 55, 137, 82, 219, 7684.0, "B-", 4300.0, 219.022136, 429.0);
    static const Species Bi_219("Bi-219", "Bi", 53, 136, 83, 219, 7700.0, "B-", 3638.0, 219.01752, 215.0);
    static const Species Po_219("Po-219", "Po", 51, 135, 84, 219, 7713.334, "B-", 2285.3395, 219.013614, 17.0);
    static const Species At_219("At-219", "At", 49, 134, 85, 219, 7720.197, "B-", 1566.6838, 219.011160587, 3.474);
    static const Species Rn_219("Rn-219", "Rn", 47, 133, 86, 219, 7723.7784, "B-", 212.3984, 219.009478683, 2.254);
    static const Species Fr_219("Fr-219", "Fr", 45, 132, 87, 219, 7721.1759, "B-", -776.9137, 219.009250664, 7.38);
    static const Species Ra_219("Ra-219", "Ra", 43, 131, 88, 219, 7714.056, "B-", -2175.7005, 219.010084715, 7.315);
    static const Species Ac_219("Ac-219", "Ac", 41, 130, 89, 219, 7700.5489, "B-", -2893.2268, 219.012420425, 55.263);
    static const Species Th_219("Th-219", "Th", 39, 129, 90, 219, 7683.7655, "B-", -4120.4434, 219.015526432, 60.611);
    static const Species Pa_219("Pa-219", "Pa", 37, 128, 91, 219, 7661.3783, "B-", -4712.7298, 219.019949909, 74.831);
    static const Species U_219("U-219", "U", 35, 127, 92, 219, 7636.2867, "B-", -6140.9976, 219.025009233, 14.319);
    static const Species Np_219("Np-219", "Np", 33, 126, 93, 219, 7604.6732, "B-", std::numeric_limits<double>::quiet_NaN(), 219.031601865, 98.732);
    static const Species Pb_220("Pb-220", "Pb", 56, 138, 82, 220, 7670.0, "B-", 3171.0, 220.025905, 429.0);
    static const Species Bi_220("Bi-220", "Bi", 54, 137, 83, 220, 7681.0, "B-", 5696.0, 220.022501, 322.0);
    static const Species Po_220("Po-220", "Po", 52, 136, 84, 220, 7703.2244, "B-", 887.7139, 220.016386, 19.0);
    static const Species At_220("At-220", "At", 50, 135, 85, 220, 7703.7033, "B-", 3763.755, 220.015433, 15.0);
    static const Species Rn_220("Rn-220", "Rn", 48, 134, 86, 220, 7717.2552, "B-", -870.3384, 220.011392443, 1.947);
    static const Species Fr_220("Fr-220", "Fr", 46, 133, 87, 220, 7709.743, "B-", 1210.2406, 220.012326789, 4.324);
    static const Species Ra_220("Ra-220", "Ra", 44, 132, 88, 220, 7711.6879, "B-", -3471.664, 220.011027542, 8.153);
    static const Species Ac_220("Ac-220", "Ac", 42, 131, 89, 220, 7692.3515, "B-", -945.7825, 220.014754527, 6.579);
    static const Species Th_220("Th-220", "Th", 40, 130, 90, 220, 7684.4964, "B-", -5588.8595, 220.015769866, 14.693);
    static const Species Pa_220("Pa-220", "Pa", 38, 129, 91, 220, 7655.5364, "B-", -2735.0, 220.021769753, 15.732);
    static const Species U_220("U-220", "U", 36, 128, 92, 220, 7640.0, "B-", -7462.0, 220.024706, 108.0);
    static const Species Np_220("Np-220", "Np", 34, 127, 93, 220, 7602.0758, "B-", std::numeric_limits<double>::quiet_NaN(), 220.03271628, 32.977);
    static const Species Bi_221("Bi-221", "Bi", 55, 138, 83, 221, 7668.0, "B-", 4426.0, 221.02598, 322.0);
    static const Species Po_221("Po-221", "Po", 53, 137, 84, 221, 7684.4814, "B-", 2991.0276, 221.021228, 21.0);
    static const Species At_221("At-221", "At", 51, 136, 85, 221, 7694.4754, "B-", 2311.375, 221.018017, 15.0);
    static const Species Rn_221("Rn-221", "Rn", 49, 135, 86, 221, 7701.3941, "B-", 1194.1032, 221.015535637, 6.134);
    static const Species Fr_221("Fr-221", "Fr", 47, 134, 87, 221, 7703.2572, "B-", 313.3741, 221.014253714, 5.245);
    static const Species Ra_221("Ra-221", "Ra", 45, 133, 88, 221, 7701.1352, "B-", -1567.1715, 221.013917293, 4.97);
    static const Species Ac_221("Ac-221", "Ac", 43, 132, 89, 221, 7690.5039, "B-", -2408.8773, 221.015599721, 61.086);
    static const Species Th_221("Th-221", "Th", 41, 131, 90, 221, 7676.064, "B-", -3435.0112, 221.018185757, 8.582);
    static const Species Pa_221("Pa-221", "Pa", 39, 130, 91, 221, 7656.9809, "B-", -4145.059, 221.021873393, 63.746);
    static const Species U_221("U-221", "U", 37, 129, 92, 221, 7634.6849, "B-", -5390.0, 221.026323297, 77.44);
    static const Species Np_221("Np-221", "Np", 35, 128, 93, 221, 7607.0, "B-", -6019.0, 221.03211, 215.0);
    static const Species Pu_221("Pu-221", "Pu", 33, 127, 94, 221, 7576.0, "B-", std::numeric_limits<double>::quiet_NaN(), 221.038572, 322.0);
    static const Species Bi_222("Bi-222", "Bi", 56, 139, 83, 222, 7648.0, "B-", 6464.0, 222.031079, 322.0);
    static const Species Po_222("Po-222", "Po", 54, 138, 84, 222, 7674.0054, "B-", 1533.2393, 222.02414, 43.0);
    static const Species At_222("At-222", "At", 52, 137, 85, 222, 7677.3878, "B-", 4581.0714, 222.022494, 17.0);
    static const Species Rn_222("Rn-222", "Rn", 50, 136, 86, 222, 7694.4991, "B-", -6.1461, 222.017576017, 2.086);
    static const Species Fr_222("Fr-222", "Fr", 48, 135, 87, 222, 7690.9474, "B-", 2057.898, 222.017582615, 8.0);
    static const Species Ra_222("Ra-222", "Ra", 46, 134, 88, 222, 7696.6931, "B-", -2301.5922, 222.015373371, 4.781);
    static const Species Ac_222("Ac-222", "Ac", 44, 133, 89, 222, 7682.8015, "B-", -581.2415, 222.017844232, 5.044);
    static const Species Th_222("Th-222", "Th", 42, 132, 90, 222, 7676.6592, "B-", -4861.322, 222.01846822, 10.966);
    static const Species Pa_222("Pa-222", "Pa", 40, 131, 91, 222, 7651.2373, "B-", -2208.4729, 222.023687064, 92.975);
    static const Species U_222("U-222", "U", 38, 130, 92, 222, 7637.7651, "B-", -7001.8072, 222.026057957, 55.817);
    static const Species Np_222("Np-222", "Np", 36, 129, 93, 222, 7602.7013, "B-", -3785.0, 222.033574706, 40.849);
    static const Species Pu_222("Pu-222", "Pu", 34, 128, 94, 222, 7582.0, "B-", std::numeric_limits<double>::quiet_NaN(), 222.037638, 322.0);
    static const Species Bi_223("Bi-223", "Bi", 57, 140, 83, 223, 7636.0, "B-", 5161.0, 223.034611, 429.0);
    static const Species Po_223("Po-223", "Po", 55, 139, 84, 223, 7655.0, "B-", 3651.0, 223.02907, 210.0);
    static const Species At_223("At-223", "At", 53, 138, 85, 223, 7668.0557, "B-", 3038.2698, 223.025151, 15.0);
    static const Species Rn_223("Rn-223", "Rn", 51, 137, 86, 223, 7678.172, "B-", 2007.4091, 223.021889283, 8.397);
    static const Species Fr_223("Fr-223", "Fr", 49, 136, 87, 223, 7683.6655, "B-", 1149.0844, 223.019734241, 2.073);
    static const Species Ra_223("Ra-223", "Ra", 47, 135, 88, 223, 7685.3101, "B-", -591.8099, 223.018500648, 2.243);
    static const Species Ac_223("Ac-223", "Ac", 45, 134, 89, 223, 7679.1479, "B-", -1560.3471, 223.019135982, 7.457);
    static const Species Th_223("Th-223", "Th", 43, 133, 90, 223, 7668.6426, "B-", -2952.2124, 223.020811083, 8.527);
    static const Species Pa_223("Pa-223", "Pa", 41, 132, 91, 223, 7651.8957, "B-", -3707.6637, 223.023980414, 81.193);
    static const Species U_223("U-223", "U", 39, 131, 92, 223, 7631.7611, "B-", -4613.3046, 223.027960754, 63.396);
    static const Species Np_223("Np-223", "Np", 37, 130, 93, 223, 7607.5654, "B-", -5462.0, 223.03291334, 88.956);
    static const Species Pu_223("Pu-223", "Pu", 35, 129, 94, 223, 7580.0, "B-", -6579.0, 223.038777, 322.0);
    static const Species Am_223("Am-223", "Am", 33, 128, 95, 223, 7547.0, "B-", std::numeric_limits<double>::quiet_NaN(), 223.04584, 322.0);
    static const Species Bi_224("Bi-224", "Bi", 58, 141, 83, 224, 7616.0, "B-", 7159.0, 224.039796, 429.0);
    static const Species Po_224("Po-224", "Po", 56, 140, 84, 224, 7644.0, "B-", 2199.0, 224.03211, 210.0);
    static const Species At_224("At-224", "At", 54, 139, 85, 224, 7650.7354, "B-", 5265.9197, 224.029749, 24.0);
    static const Species Rn_224("Rn-224", "Rn", 52, 138, 86, 224, 7670.7514, "B-", 696.484, 224.024095803, 10.536);
    static const Species Fr_224("Fr-224", "Fr", 50, 137, 87, 224, 7670.368, "B-", 2922.7819, 224.023348096, 12.0);
    static const Species Ra_224("Ra-224", "Ra", 48, 136, 88, 224, 7679.9236, "B-", -1408.3152, 224.020210361, 1.944);
    static const Species Ac_224("Ac-224", "Ac", 46, 135, 89, 224, 7670.1438, "B-", 238.5672, 224.021722249, 4.389);
    static const Species Th_224("Th-224", "Th", 44, 134, 90, 224, 7667.7162, "B-", -3866.7705, 224.021466137, 10.31);
    static const Species Pa_224("Pa-224", "Pa", 42, 133, 91, 224, 7646.9612, "B-", -1880.3393, 224.025617286, 8.145);
    static const Species U_224("U-224", "U", 40, 132, 92, 224, 7635.0742, "B-", -6289.5572, 224.027635913, 16.383);
    static const Species Np_224("Np-224", "Np", 38, 131, 93, 224, 7603.5032, "B-", -3248.0, 224.03438803, 31.052);
    static const Species Pu_224("Pu-224", "Pu", 36, 130, 94, 224, 7586.0, "B-", -7980.0, 224.037875, 322.0);
    static const Species Am_224("Am-224", "Am", 34, 129, 95, 224, 7546.0, "B-", std::numeric_limits<double>::quiet_NaN(), 224.046442, 429.0);
    static const Species Po_225("Po-225", "Po", 57, 141, 84, 225, 7626.0, "B-", 4280.0, 225.037123, 322.0);
    static const Species At_225("At-225", "At", 55, 140, 85, 225, 7641.0, "B-", 3765.0, 225.032528, 322.0);
    static const Species Rn_225("Rn-225", "Rn", 53, 139, 86, 225, 7654.3581, "B-", 2713.5412, 225.028485572, 11.958);
    static const Species Fr_225("Fr-225", "Fr", 51, 138, 87, 225, 7662.9412, "B-", 1827.5584, 225.025572466, 12.847);
    static const Species Ra_225("Ra-225", "Ra", 49, 137, 88, 225, 7667.5866, "B-", 355.7386, 225.023610502, 2.786);
    static const Species Ac_225("Ac-225", "Ac", 47, 136, 89, 225, 7665.6906, "B-", -672.8878, 225.023228601, 5.107);
    static const Species Th_225("Th-225", "Th", 45, 135, 90, 225, 7659.2229, "B-", -2046.4473, 225.023950975, 5.467);
    static const Species Pa_225("Pa-225", "Pa", 43, 134, 91, 225, 7646.6504, "B-", -3015.361, 225.026147927, 87.887);
    static const Species U_225("U-225", "U", 41, 133, 92, 225, 7629.7717, "B-", -4246.0969, 225.02938505, 10.664);
    static const Species Np_225("Np-225", "Np", 39, 132, 93, 225, 7607.4231, "B-", -4682.0, 225.033943422, 98.355);
    static const Species Pu_225("Pu-225", "Pu", 37, 131, 94, 225, 7583.0, "B-", -6090.0, 225.03897, 322.0);
    static const Species Am_225("Am-225", "Am", 35, 130, 95, 225, 7553.0, "B-", std::numeric_limits<double>::quiet_NaN(), 225.045508, 429.0);
    static const Species Po_226("Po-226", "Po", 58, 142, 84, 226, 7614.0, "B-", 2889.0, 226.04031, 430.0);
    static const Species At_226("At-226", "At", 56, 141, 85, 226, 7624.0, "B-", 5913.0, 226.037209, 322.0);
    static const Species Rn_226("Rn-226", "Rn", 54, 140, 86, 226, 7646.4108, "B-", 1226.6542, 226.03086138, 11.247);
    static const Species Fr_226("Fr-226", "Fr", 52, 139, 87, 226, 7648.3768, "B-", 3852.9638, 226.029544512, 6.688);
    static const Species Ra_226("Ra-226", "Ra", 50, 138, 88, 226, 7661.9636, "B-", -641.6252, 226.025408186, 2.068);
    static const Species Ac_226("Ac-226", "Ac", 48, 137, 89, 226, 7655.6628, "B-", 1111.5517, 226.026096999, 3.327);
    static const Species Th_226("Th-226", "Th", 46, 136, 90, 226, 7657.1195, "B-", -2835.9504, 226.024903699, 4.81);
    static const Species Pa_226("Pa-226", "Pa", 44, 135, 91, 226, 7641.1093, "B-", -1295.1978, 226.027948217, 12.037);
    static const Species U_226("U-226", "U", 42, 134, 92, 226, 7631.9166, "B-", -5488.0792, 226.029338669, 11.884);
    static const Species Np_226("Np-226", "Np", 40, 133, 93, 226, 7604.1714, "B-", -2813.0, 226.035230364, 109.568);
    static const Species Pu_226("Pu-226", "Pu", 38, 132, 94, 226, 7588.0, "B-", -7340.0, 226.03825, 215.0);
    static const Species Am_226("Am-226", "Am", 36, 131, 95, 226, 7552.0, "B-", std::numeric_limits<double>::quiet_NaN(), 226.04613, 322.0);
    static const Species Po_227("Po-227", "Po", 59, 143, 84, 227, 7596.0, "B-", 4850.0, 227.04539, 430.0);
    static const Species At_227("At-227", "At", 57, 142, 85, 227, 7613.0, "B-", 4544.0, 227.040183, 322.0);
    static const Species Rn_227("Rn-227", "Rn", 55, 141, 86, 227, 7630.0508, "B-", 3203.3894, 227.035304393, 15.127);
    static const Species Fr_227("Fr-227", "Fr", 53, 140, 87, 227, 7640.7162, "B-", 2504.9813, 227.031865413, 6.332);
    static const Species Ra_227("Ra-227", "Ra", 51, 139, 88, 227, 7648.3048, "B-", 1327.9489, 227.029176205, 2.089);
    static const Species Ac_227("Ac-227", "Ac", 49, 138, 89, 227, 7650.7084, "B-", 44.7559, 227.027750594, 2.068);
    static const Species Th_227("Th-227", "Th", 47, 137, 90, 227, 7647.4591, "B-", -1025.6117, 227.027702546, 2.241);
    static const Species Pa_227("Pa-227", "Pa", 45, 136, 91, 227, 7639.4945, "B-", -2214.6629, 227.028803586, 7.797);
    static const Species U_227("U-227", "U", 43, 135, 92, 227, 7626.2918, "B-", -3533.9848, 227.031181124, 9.136);
    static const Species Np_227("Np-227", "Np", 41, 134, 93, 227, 7607.2771, "B-", -4191.0, 227.034975012, 82.651);
    static const Species Pu_227("Pu-227", "Pu", 39, 133, 94, 227, 7585.0, "B-", -5410.0, 227.039474, 107.0);
    static const Species Am_227("Am-227", "Am", 37, 132, 95, 227, 7558.0, "B-", std::numeric_limits<double>::quiet_NaN(), 227.045282, 215.0);
    static const Species At_228("At-228", "At", 58, 143, 85, 228, 7596.0, "B-", 6637.0, 228.04496, 429.0);
    static const Species Rn_228("Rn-228", "Rn", 56, 142, 86, 228, 7621.6457, "B-", 1859.2451, 228.037835415, 18.977);
    static const Species Fr_228("Fr-228", "Fr", 54, 141, 87, 228, 7626.3689, "B-", 4444.027, 228.035839433, 7.226);
    static const Species Ra_228("Ra-228", "Ra", 52, 140, 88, 228, 7642.4289, "B-", 45.5402, 228.031068574, 2.141);
    static const Species Ac_228("Ac-228", "Ac", 50, 139, 89, 228, 7639.1973, "B-", 2123.7545, 228.031019685, 2.247);
    static const Species Th_228("Th-228", "Th", 48, 138, 90, 228, 7645.0807, "B-", -2152.6993, 228.028739741, 1.938);
    static const Species Pa_228("Pa-228", "Pa", 46, 137, 91, 228, 7632.2076, "B-", -296.402, 228.031050758, 4.659);
    static const Species U_228("U-228", "U", 44, 136, 92, 228, 7627.4763, "B-", -4605.0, 228.031368959, 14.465);
    static const Species Np_228("Np-228", "Np", 42, 135, 93, 228, 7604.0, "B-", -2283.0, 228.036313, 108.0);
    static const Species Pu_228("Pu-228", "Pu", 40, 134, 94, 228, 7590.4039, "B-", -6742.0, 228.038763325, 25.069);
    static const Species Am_228("Am-228", "Am", 38, 133, 95, 228, 7557.0, "B-", std::numeric_limits<double>::quiet_NaN(), 228.046001, 215.0);
    static const Species At_229("At-229", "At", 59, 144, 85, 229, 7585.0, "B-", 5527.0, 229.048191, 429.0);
    static const Species Rn_229("Rn-229", "Rn", 57, 143, 86, 229, 7605.6227, "B-", 3694.1465, 229.042257272, 14.0);
    static const Species Fr_229("Fr-229", "Fr", 55, 142, 87, 229, 7618.338, "B-", 3106.2907, 229.038291443, 5.368);
    static const Species Ra_229("Ra-229", "Ra", 53, 141, 88, 229, 7628.4862, "B-", 1872.0266, 229.034956703, 16.576);
    static const Species Ac_229("Ac-229", "Ac", 51, 140, 89, 229, 7633.2446, "B-", 1104.4191, 229.032947, 13.0);
    static const Species Th_229("Th-229", "Th", 49, 139, 90, 229, 7634.651, "B-", -311.331, 229.031761357, 2.581);
    static const Species Pa_229("Pa-229", "Pa", 47, 138, 91, 229, 7629.8752, "B-", -1313.7716, 229.032095585, 3.521);
    static const Species U_229("U-229", "U", 45, 137, 92, 229, 7620.7218, "B-", -2590.7577, 229.033505976, 6.374);
    static const Species Np_229("Np-229", "Np", 43, 136, 93, 229, 7605.9921, "B-", -3593.5462, 229.036287269, 108.618);
    static const Species Pu_229("Pu-229", "Pu", 41, 135, 94, 229, 7586.8834, "B-", -4785.4899, 229.040145099, 65.092);
    static const Species Am_229("Am-229", "Am", 39, 134, 95, 229, 7562.5697, "B-", std::numeric_limits<double>::quiet_NaN(), 229.045282534, 114.169);
    static const Species Rn_230("Rn-230", "Rn", 58, 144, 86, 230, 7595.0, "B-", 2683.0, 230.045271, 215.0);
    static const Species Fr_230("Fr-230", "Fr", 56, 143, 87, 230, 7603.7052, "B-", 4970.4627, 230.042390787, 7.022);
    static const Species Ra_230("Ra-230", "Ra", 54, 142, 88, 230, 7621.9144, "B-", 677.9196, 230.037054776, 11.053);
    static const Species Ac_230("Ac-230", "Ac", 52, 141, 89, 230, 7621.4604, "B-", 2975.8745, 230.036327, 17.0);
    static const Species Th_230("Th-230", "Th", 50, 140, 90, 230, 7630.9974, "B-", -1311.0313, 230.033132267, 1.297);
    static const Species Pa_230("Pa-230", "Pa", 48, 139, 91, 230, 7621.8958, "B-", 558.5262, 230.034539717, 3.261);
    static const Species U_230("U-230", "U", 46, 138, 92, 230, 7620.9227, "B-", -3621.5986, 230.033940114, 4.841);
    static const Species Np_230("Np-230", "Np", 44, 137, 93, 230, 7601.7751, "B-", -1695.5543, 230.03782806, 59.051);
    static const Species Pu_230("Pu-230", "Pu", 42, 136, 94, 230, 7591.0016, "B-", -5940.0, 230.039648313, 15.514);
    static const Species Am_230("Am-230", "Am", 40, 135, 95, 230, 7562.0, "B-", std::numeric_limits<double>::quiet_NaN(), 230.046025, 153.0);
    static const Species Rn_231("Rn-231", "Rn", 59, 145, 86, 231, 7579.0, "B-", 4469.0, 231.049973, 322.0);
    static const Species Fr_231("Fr-231", "Fr", 57, 144, 87, 231, 7594.5009, "B-", 3864.0868, 231.045175353, 8.3);
    static const Species Ra_231("Ra-231", "Ra", 55, 143, 88, 231, 7607.8418, "B-", 2453.6351, 231.041027085, 12.206);
    static const Species Ac_231("Ac-231", "Ac", 53, 142, 89, 231, 7615.0768, "B-", 1947.0425, 231.038393, 14.0);
    static const Species Th_231("Th-231", "Th", 51, 141, 90, 231, 7620.1188, "B-", 391.4727, 231.036302764, 1.306);
    static const Species Pa_231("Pa-231", "Pa", 49, 140, 91, 231, 7618.4267, "B-", -381.6138, 231.0358825, 1.901);
    static const Species U_231("U-231", "U", 47, 139, 92, 231, 7613.3879, "B-", -1817.7347, 231.03629218, 2.866);
    static const Species Np_231("Np-231", "Np", 45, 138, 93, 231, 7602.1321, "B-", -2684.8905, 231.038243598, 54.916);
    static const Species Pu_231("Pu-231", "Pu", 43, 137, 94, 231, 7587.1224, "B-", -4101.0, 231.041125946, 23.683);
    static const Species Am_231("Am-231", "Am", 41, 136, 95, 231, 7566.0, "B-", -4860.0, 231.045529, 322.0);
    static const Species Cm_231("Cm-231", "Cm", 39, 135, 96, 231, 7542.0, "B-", std::numeric_limits<double>::quiet_NaN(), 231.050746, 322.0);
    static const Species Fr_232("Fr-232", "Fr", 58, 145, 87, 232, 7579.3481, "B-", 5575.8791, 232.049461219, 15.0);
    static const Species Ra_232("Ra-232", "Ra", 56, 144, 88, 232, 7600.0099, "B-", 1342.5322, 232.043475267, 9.823);
    static const Species Ac_232("Ac-232", "Ac", 54, 143, 89, 232, 7602.4245, "B-", 3707.7131, 232.042034, 14.0);
    static const Species Th_232("Th-232", "Th", 52, 142, 90, 232, 7615.0338, "B-", -499.8388, 232.038053606, 1.525);
    static const Species Pa_232("Pa-232", "Pa", 50, 141, 91, 232, 7609.5072, "B-", 1337.1034, 232.038590205, 8.206);
    static const Species U_232("U-232", "U", 48, 140, 92, 232, 7611.8984, "B-", -2750.0, 232.037154765, 1.941);
    static const Species Np_232("Np-232", "Np", 46, 139, 93, 232, 7597.0, "B-", -1001.0, 232.040107, 107.0);
    static const Species Pu_232("Pu-232", "Pu", 44, 138, 94, 232, 7588.9839, "B-", -5059.0, 232.041182133, 18.126);
    static const Species Am_232("Am-232", "Am", 42, 137, 95, 232, 7564.0, "B-", -2913.0, 232.046613, 322.0);
    static const Species Cm_232("Cm-232", "Cm", 40, 136, 96, 232, 7548.0, "B-", std::numeric_limits<double>::quiet_NaN(), 232.04974, 216.0);
    static const Species Fr_233("Fr-233", "Fr", 59, 146, 87, 233, 7569.2398, "B-", 4585.9906, 233.052517833, 21.0);
    static const Species Ra_233("Ra-233", "Ra", 57, 145, 88, 233, 7585.5644, "B-", 3026.0244, 233.04759457, 9.235);
    static const Species Ac_233("Ac-233", "Ac", 55, 144, 89, 233, 7595.1939, "B-", 2576.395, 233.044346, 14.0);
    static const Species Th_233("Th-233", "Th", 53, 143, 90, 233, 7602.8937, "B-", 1242.232, 233.041580126, 1.528);
    static const Species Pa_233("Pa-233", "Pa", 51, 142, 91, 233, 7604.8675, "B-", 570.2993, 233.040246535, 1.433);
    static const Species U_233("U-233", "U", 49, 141, 92, 233, 7603.9574, "B-", -1029.4197, 233.039634294, 2.42);
    static const Species Np_233("Np-233", "Np", 47, 140, 93, 233, 7596.1816, "B-", -2103.3047, 233.040739421, 54.729);
    static const Species Pu_233("Pu-233", "Pu", 45, 139, 94, 233, 7583.7968, "B-", -3233.0, 233.042997411, 58.162);
    static const Species Am_233("Am-233", "Am", 43, 138, 95, 233, 7567.0, "B-", -4008.0, 233.046468, 123.0);
    static const Species Cm_233("Cm-233", "Cm", 41, 137, 96, 233, 7546.002, "B-", -5478.0, 233.050771485, 87.059);
    static const Species Bk_233("Bk-233", "Bk", 39, 136, 97, 233, 7519.0, "B-", std::numeric_limits<double>::quiet_NaN(), 233.056652, 250.0);
    static const Species Ra_234("Ra-234", "Ra", 58, 146, 88, 234, 7576.5439, "B-", 2089.4348, 234.0503821, 9.0);
    static const Species Ac_234("Ac-234", "Ac", 56, 145, 89, 234, 7582.1297, "B-", 4228.2364, 234.048139, 15.0);
    static const Species Th_234("Th-234", "Th", 54, 144, 90, 234, 7596.8557, "B-", 274.0882, 234.043599801, 2.779);
    static const Species Pa_234("Pa-234", "Pa", 52, 143, 91, 234, 7594.6837, "B-", 2193.9105, 234.043305555, 4.395);
    static const Species U_234("U-234", "U", 50, 142, 92, 234, 7600.716, "B-", -1809.8462, 234.040950296, 1.212);
    static const Species Np_234("Np-234", "Np", 48, 141, 93, 234, 7589.6382, "B-", -395.1807, 234.042893245, 9.014);
    static const Species Pu_234("Pu-234", "Pu", 46, 140, 94, 234, 7584.6061, "B-", -4112.0, 234.043317489, 7.298);
    static const Species Am_234("Am-234", "Am", 44, 139, 95, 234, 7564.0, "B-", -2261.0, 234.047731, 172.0);
    static const Species Cm_234("Cm-234", "Cm", 42, 138, 96, 234, 7550.6868, "B-", -6673.0, 234.050158568, 18.333);
    static const Species Bk_234("Bk-234", "Bk", 40, 137, 97, 234, 7519.0, "B-", std::numeric_limits<double>::quiet_NaN(), 234.057322, 164.0);
    static const Species Ra_235("Ra-235", "Ra", 59, 147, 88, 235, 7561.0, "B-", 3773.0, 235.05489, 322.0);
    static const Species Ac_235("Ac-235", "Ac", 57, 146, 89, 235, 7573.5051, "B-", 3339.4064, 235.05084, 15.0);
    static const Species Th_235("Th-235", "Th", 55, 145, 90, 235, 7584.3862, "B-", 1728.8531, 235.047255, 14.0);
    static const Species Pa_235("Pa-235", "Pa", 53, 144, 91, 235, 7588.4139, "B-", 1370.1184, 235.045399, 15.0);
    static const Species U_235("U-235", "U", 51, 143, 92, 235, 7590.9151, "B-", -124.2619, 235.043928117, 1.198);
    static const Species Np_235("Np-235", "Np", 49, 142, 93, 235, 7587.0571, "B-", -1139.3021, 235.044061518, 1.49);
    static const Species Pu_235("Pu-235", "Pu", 47, 141, 94, 235, 7578.8799, "B-", -2442.2558, 235.045284609, 22.03);
    static const Species Am_235("Am-235", "Am", 45, 140, 95, 235, 7565.1582, "B-", -3389.0, 235.047906478, 56.661);
    static const Species Cm_235("Cm-235", "Cm", 43, 139, 96, 235, 7547.0, "B-", -4757.0, 235.051545, 110.0);
    static const Species Bk_235("Bk-235", "Bk", 41, 138, 97, 235, 7524.0, "B-", std::numeric_limits<double>::quiet_NaN(), 235.056651, 430.0);
    static const Species Ac_236("Ac-236", "Ac", 58, 147, 89, 236, 7559.2423, "B-", 4965.7951, 236.054988, 41.0);
    static const Species Th_236("Th-236", "Th", 56, 146, 90, 236, 7576.9688, "B-", 921.2477, 236.049657, 15.0);
    static const Species Pa_236("Pa-236", "Pa", 54, 145, 91, 236, 7577.5573, "B-", 2889.373, 236.048668, 15.0);
    static const Species U_236("U-236", "U", 52, 144, 92, 236, 7586.4854, "B-", -933.5116, 236.04556613, 1.193);
    static const Species Np_236("Np-236", "Np", 50, 143, 93, 236, 7579.2148, "B-", 476.5854, 236.046568296, 54.129);
    static const Species Pu_236("Pu-236", "Pu", 48, 142, 94, 236, 7577.9192, "B-", -3139.0, 236.046056661, 1.942);
    static const Species Am_236("Am-236", "Am", 46, 141, 95, 236, 7561.0, "B-", -1812.0, 236.049427, 127.0);
    static const Species Cm_236("Cm-236", "Cm", 44, 140, 96, 236, 7550.309, "B-", -5689.0, 236.051372112, 18.931);
    static const Species Bk_236("Bk-236", "Bk", 42, 139, 97, 236, 7523.0, "B-", std::numeric_limits<double>::quiet_NaN(), 236.057479, 387.0);
    static const Species Ac_237("Ac-237", "Ac", 59, 148, 89, 237, 7550.0, "B-", 4065.0, 237.057993, 429.0);
    static const Species Th_237("Th-237", "Th", 57, 147, 90, 237, 7563.4433, "B-", 2427.4736, 237.053629, 17.0);
    static const Species Pa_237("Pa-237", "Pa", 55, 146, 91, 237, 7570.3847, "B-", 2137.4905, 237.051023, 14.0);
    static const Species U_237("U-237", "U", 53, 145, 92, 237, 7576.1026, "B-", 518.5338, 237.048728309, 1.29);
    static const Species Np_237("Np-237", "Np", 51, 144, 93, 237, 7574.9895, "B-", -220.063, 237.04817164, 1.201);
    static const Species Pu_237("Pu-237", "Pu", 49, 143, 94, 237, 7570.7599, "B-", -1478.0, 237.048407888, 1.821);
    static const Species Am_237("Am-237", "Am", 47, 142, 95, 237, 7561.0, "B-", -2677.0, 237.049995, 64.0);
    static const Species Cm_237("Cm-237", "Cm", 45, 141, 96, 237, 7546.6241, "B-", -3963.0, 237.052868988, 79.87);
    static const Species Bk_237("Bk-237", "Bk", 43, 140, 97, 237, 7527.0, "B-", -4728.0, 237.057123, 247.0);
    static const Species Cf_237("Cf-237", "Cf", 41, 139, 98, 237, 7503.3507, "B-", std::numeric_limits<double>::quiet_NaN(), 237.062199272, 104.506);
    static const Species Th_238("Th-238", "Th", 58, 148, 90, 238, 7555.0, "B-", 1631.0, 238.056388, 304.0);
    static const Species Pa_238("Pa-238", "Pa", 56, 147, 91, 238, 7558.3449, "B-", 3586.3111, 238.054637, 17.0);
    static const Species U_238("U-238", "U", 54, 146, 92, 238, 7570.1262, "B-", -146.8652, 238.050786936, 1.601);
    static const Species Np_238("Np-238", "Np", 52, 145, 93, 238, 7566.222, "B-", 1291.4491, 238.050944603, 1.22);
    static const Species Pu_238("Pu-238", "Pu", 50, 144, 94, 238, 7568.3611, "B-", -2258.2731, 238.049558175, 1.221);
    static const Species Am_238("Am-238", "Am", 48, 143, 95, 238, 7555.5853, "B-", -1023.7818, 238.051982531, 63.243);
    static const Species Cm_238("Cm-238", "Cm", 46, 142, 96, 238, 7547.9966, "B-", -4771.0, 238.053081606, 13.133);
    static const Species Bk_238("Bk-238", "Bk", 44, 141, 97, 238, 7525.0, "B-", -3061.0, 238.058204, 275.0);
    static const Species Cf_238("Cf-238", "Cf", 42, 140, 98, 238, 7509.0, "B-", std::numeric_limits<double>::quiet_NaN(), 238.06149, 320.0);
    static const Species Th_239("Th-239", "Th", 59, 149, 90, 239, 7540.0, "B-", 3162.0, 239.060655, 429.0);
    static const Species Pa_239("Pa-239", "Pa", 57, 148, 91, 239, 7550.0, "B-", 2765.0, 239.05726, 210.0);
    static const Species U_239("U-239", "U", 55, 147, 92, 239, 7558.5624, "B-", 1261.6634, 239.054291989, 1.612);
    static const Species Np_239("Np-239", "Np", 53, 146, 93, 239, 7560.568, "B-", 722.7849, 239.052937538, 1.406);
    static const Species Pu_239("Pu-239", "Pu", 51, 145, 94, 239, 7560.3187, "B-", -802.1402, 239.052161596, 1.194);
    static const Species Am_239("Am-239", "Am", 49, 144, 95, 239, 7553.6891, "B-", -1756.6021, 239.053022729, 2.127);
    static const Species Cm_239("Cm-239", "Cm", 47, 143, 96, 239, 7543.0659, "B-", -3103.0, 239.054908519, 161.107);
    static const Species Bk_239("Bk-239", "Bk", 45, 142, 97, 239, 7527.0, "B-", -3952.0, 239.058239, 222.0);
    static const Species Cf_239("Cf-239", "Cf", 43, 141, 98, 239, 7507.0, "B-", -5429.0, 239.062482, 129.0);
    static const Species Es_239("Es-239", "Es", 41, 140, 99, 239, 7481.0, "B-", std::numeric_limits<double>::quiet_NaN(), 239.06831, 322.0);
    static const Species Pa_240("Pa-240", "Pa", 58, 149, 91, 240, 7537.0, "B-", 4295.0, 240.061203, 215.0);
    static const Species U_240("U-240", "U", 56, 148, 92, 240, 7551.7705, "B-", 399.2685, 240.056592411, 2.74);
    static const Species Np_240("Np-240", "Np", 54, 147, 93, 240, 7550.1743, "B-", 2190.9095, 240.056163778, 18.284);
    static const Species Pu_240("Pu-240", "Pu", 52, 146, 94, 240, 7556.0433, "B-", -1384.7902, 240.05381174, 1.186);
    static const Species Am_240("Am-240", "Am", 50, 145, 95, 240, 7547.0136, "B-", -214.1127, 240.055298374, 14.849);
    static const Species Cm_240("Cm-240", "Cm", 48, 144, 96, 240, 7542.8617, "B-", -3940.0, 240.055528233, 2.045);
    static const Species Bk_240("Bk-240", "Bk", 46, 143, 97, 240, 7523.0, "B-", -2324.0, 240.059758, 161.0);
    static const Species Cf_240("Cf-240", "Cf", 44, 142, 98, 240, 7510.24, "B-", -6237.0, 240.062253447, 19.36);
    static const Species Es_240("Es-240", "Es", 42, 141, 99, 240, 7481.0, "B-", std::numeric_limits<double>::quiet_NaN(), 240.068949, 393.0);
    static const Species Pa_241("Pa-241", "Pa", 59, 150, 91, 241, 7528.0, "B-", 3543.0, 241.064134, 322.0);
    static const Species U_241("U-241", "U", 57, 149, 92, 241, 7539.0, "B-", 1882.0, 241.06033, 210.0);
    static const Species Np_241("Np-241", "Np", 55, 148, 93, 241, 7544.0426, "B-", 1360.0, 241.058309671, 107.36);
    static const Species Pu_241("Pu-241", "Pu", 53, 147, 94, 241, 7546.4395, "B-", 20.7799, 241.056849651, 1.186);
    static const Species Am_241("Am-241", "Am", 51, 146, 95, 241, 7543.2795, "B-", -767.4346, 241.056827343, 1.195);
    static const Species Cm_241("Cm-241", "Cm", 49, 145, 96, 241, 7536.8488, "B-", -2279.0, 241.057651218, 1.725);
    static const Species Bk_241("Bk-241", "Bk", 47, 144, 97, 241, 7524.0, "B-", -3346.0, 241.060098, 178.0);
    static const Species Cf_241("Cf-241", "Cf", 45, 143, 98, 241, 7507.0, "B-", -4567.0, 241.06369, 180.0);
    static const Species Es_241("Es-241", "Es", 43, 142, 99, 241, 7485.0, "B-", -5327.0, 241.068592, 248.0);
    static const Species Fm_241("Fm-241", "Fm", 41, 141, 100, 241, 7459.0, "B-", std::numeric_limits<double>::quiet_NaN(), 241.074311, 322.0);
    static const Species U_242("U-242", "U", 58, 150, 92, 242, 7532.0, "B-", 1203.0, 242.062931, 215.0);
    static const Species Np_242("Np-242", "Np", 56, 149, 93, 242, 7533.4042, "B-", 2700.0, 242.061639548, 214.712);
    static const Species Pu_242("Pu-242", "Pu", 54, 148, 94, 242, 7541.3284, "B-", -751.1373, 242.058740979, 1.336);
    static const Species Am_242("Am-242", "Am", 52, 147, 95, 242, 7534.9917, "B-", 664.3145, 242.059547358, 1.199);
    static const Species Cm_242("Cm-242", "Cm", 50, 146, 96, 242, 7534.504, "B-", -2948.0, 242.058834187, 1.224);
    static const Species Bk_242("Bk-242", "Bk", 48, 145, 97, 242, 7519.0, "B-", -1635.0, 242.061999, 144.0);
    static const Species Cf_242("Cf-242", "Cf", 46, 144, 98, 242, 7509.0991, "B-", -5414.0, 242.063754544, 13.84);
    static const Species Es_242("Es-242", "Es", 44, 143, 99, 242, 7483.0, "B-", -3598.0, 242.069567, 276.0);
    static const Species Fm_242("Fm-242", "Fm", 42, 142, 100, 242, 7465.0, "B-", std::numeric_limits<double>::quiet_NaN(), 242.07343, 430.0);
    static const Species U_243("U-243", "U", 59, 151, 92, 243, 7518.0, "B-", 2674.0, 243.067075, 322.0);
    static const Species Np_243("Np-243", "Np", 57, 150, 93, 243, 7526.0, "B-", 2051.0, 243.064204, 34.0);
    static const Species Pu_243("Pu-243", "Pu", 55, 149, 94, 243, 7531.0087, "B-", 579.5559, 243.062002068, 2.728);
    static const Species Am_243("Am-243", "Am", 53, 148, 95, 243, 7530.1742, "B-", -6.9302, 243.061379889, 1.49);
    static const Species Cm_243("Cm-243", "Cm", 51, 147, 96, 243, 7526.9261, "B-", -1507.6936, 243.061387329, 1.605);
    static const Species Bk_243("Bk-243", "Bk", 49, 146, 97, 243, 7517.5021, "B-", -2300.0, 243.063005905, 4.856);
    static const Species Cf_243("Cf-243", "Cf", 47, 145, 98, 243, 7505.0, "B-", -3757.0, 243.065475, 194.0);
    static const Species Es_243("Es-243", "Es", 45, 144, 99, 243, 7486.0, "B-", -4569.0, 243.069508, 222.0);
    static const Species Fm_243("Fm-243", "Fm", 43, 143, 100, 243, 7464.0, "B-", std::numeric_limits<double>::quiet_NaN(), 243.074414, 140.0);
    static const Species Np_244("Np-244", "Np", 58, 151, 93, 244, 7514.0, "B-", 3434.0, 244.067891, 107.0);
    static const Species Pu_244("Pu-244", "Pu", 56, 150, 94, 244, 7524.8154, "B-", -73.1143, 244.064204401, 2.518);
    static const Species Am_244("Am-244", "Am", 54, 149, 95, 244, 7521.3095, "B-", 1427.3, 244.064282892, 1.6);
    static const Species Cm_244("Cm-244", "Cm", 52, 148, 96, 244, 7523.9527, "B-", -2261.9902, 244.062750622, 1.187);
    static const Species Bk_244("Bk-244", "Bk", 50, 147, 97, 244, 7511.4759, "B-", -764.2709, 244.065178969, 15.457);
    static const Species Cf_244("Cf-244", "Cf", 48, 146, 98, 244, 7505.1373, "B-", -4547.0, 244.065999447, 2.809);
    static const Species Es_244("Es-244", "Es", 46, 145, 99, 244, 7483.0, "B-", -2938.0, 244.070881, 195.0);
    static const Species Fm_244("Fm-244", "Fm", 44, 144, 100, 244, 7468.0, "B-", -6634.0, 244.074036, 216.0);
    static const Species Md_244("Md-244", "Md", 42, 143, 101, 244, 7438.0, "B-", std::numeric_limits<double>::quiet_NaN(), 244.081157, 402.0);
    static const Species Np_245("Np-245", "Np", 59, 152, 93, 245, 7506.0, "B-", 2672.0, 245.070693, 215.0);
    static const Species Pu_245("Pu-245", "Pu", 57, 151, 94, 245, 7513.2822, "B-", 1277.7559, 245.067824554, 14.621);
    static const Species Am_245("Am-245", "Am", 55, 150, 95, 245, 7515.3043, "B-", 895.8929, 245.066452827, 2.024);
    static const Species Cm_245("Cm-245", "Cm", 53, 149, 96, 245, 7515.7677, "B-", -809.2519, 245.065491047, 1.233);
    static const Species Bk_245("Bk-245", "Bk", 51, 148, 97, 245, 7509.2714, "B-", -1571.3755, 245.066359814, 1.923);
    static const Species Cf_245("Cf-245", "Cf", 49, 147, 98, 245, 7499.6644, "B-", -2930.0, 245.068046755, 2.606);
    static const Species Es_245("Es-245", "Es", 47, 146, 99, 245, 7485.0, "B-", -3877.0, 245.071192, 178.0);
    static const Species Fm_245("Fm-245", "Fm", 45, 145, 100, 245, 7465.0, "B-", -5133.0, 245.075354, 210.0);
    static const Species Md_245("Md-245", "Md", 43, 144, 101, 245, 7441.0, "B-", std::numeric_limits<double>::quiet_NaN(), 245.080864, 279.0);
    static const Species Pu_246("Pu-246", "Pu", 58, 152, 94, 246, 7506.5401, "B-", 401.0, 246.070204172, 16.087);
    static const Species Am_246("Am-246", "Am", 56, 151, 95, 246, 7505.0, "B-", 2377.0, 246.069774, 19.0);
    static const Species Cm_246("Cm-246", "Cm", 54, 150, 96, 246, 7511.4716, "B-", -1350.0, 246.067222016, 1.637);
    static const Species Bk_246("Bk-246", "Bk", 52, 149, 97, 246, 7502.8035, "B-", -123.3159, 246.0686713, 64.433);
    static const Species Cf_246("Cf-246", "Cf", 50, 148, 98, 246, 7499.122, "B-", -3728.5741, 246.068803685, 1.625);
    static const Species Es_246("Es-246", "Es", 48, 147, 99, 246, 7480.7849, "B-", -2372.3848, 246.072806474, 96.538);
    static const Species Fm_246("Fm-246", "Fm", 46, 146, 100, 246, 7467.9608, "B-", -5924.0, 246.075353334, 14.675);
    static const Species Md_246("Md-246", "Md", 44, 145, 101, 246, 7441.0, "B-", std::numeric_limits<double>::quiet_NaN(), 246.081713, 279.0);
    static const Species Pu_247("Pu-247", "Pu", 59, 153, 94, 247, 7493.0, "B-", 2057.0, 247.0743, 215.0);
    static const Species Am_247("Am-247", "Am", 57, 152, 95, 247, 7499.0, "B-", 1620.0, 247.072092, 107.0);
    static const Species Cm_247("Cm-247", "Cm", 55, 151, 96, 247, 7501.9318, "B-", 43.5841, 247.070352678, 4.076);
    static const Species Bk_247("Bk-247", "Bk", 53, 150, 97, 247, 7498.9408, "B-", -619.8711, 247.070305889, 5.57);
    static const Species Cf_247("Cf-247", "Cf", 51, 149, 98, 247, 7493.2638, "B-", -2469.0006, 247.070971348, 15.38);
    static const Species Es_247("Es-247", "Es", 49, 148, 99, 247, 7480.1005, "B-", -3094.0, 247.073621929, 20.87);
    static const Species Fm_247("Fm-247", "Fm", 47, 147, 100, 247, 7464.0, "B-", -4263.0, 247.076944, 194.0);
    static const Species Md_247("Md-247", "Md", 45, 146, 101, 247, 7444.0, "B-", std::numeric_limits<double>::quiet_NaN(), 247.08152, 223.0);
    static const Species Am_248("Am-248", "Am", 58, 153, 95, 248, 7487.0, "B-", 3170.0, 248.075752, 215.0);
    static const Species Cm_248("Cm-248", "Cm", 56, 152, 96, 248, 7496.7291, "B-", -738.3049, 248.072349086, 2.531);
    static const Species Bk_248("Bk-248", "Bk", 54, 151, 97, 248, 7490.5975, "B-", 893.1015, 248.073141689, 53.739);
    static const Species Cf_248("Cf-248", "Cf", 52, 150, 98, 248, 7491.044, "B-", -3061.0, 248.072182905, 5.497);
    static const Species Es_248("Es-248", "Es", 50, 149, 99, 248, 7476.0, "B-", -1599.0, 248.075469, 56.0);
    static const Species Fm_248("Fm-248", "Fm", 48, 148, 100, 248, 7465.9451, "B-", -5050.0, 248.077185451, 9.122);
    static const Species Md_248("Md-248", "Md", 46, 147, 101, 248, 7442.0, "B-", -3741.0, 248.082607, 198.0);
    static const Species No_248("No-248", "No", 44, 146, 102, 248, 7424.0, "B-", std::numeric_limits<double>::quiet_NaN(), 248.086623, 241.0);
    static const Species Am_249("Am-249", "Am", 59, 154, 95, 249, 7479.0, "B-", 2353.0, 249.07848, 320.0);
    static const Species Cm_249("Cm-249", "Cm", 57, 153, 96, 249, 7485.551, "B-", 904.363, 249.075953992, 2.545);
    static const Species Bk_249("Bk-249", "Bk", 55, 152, 97, 249, 7486.041, "B-", 123.6, 249.074983118, 1.339);
    static const Species Cf_249("Cf-249", "Cf", 53, 151, 98, 249, 7483.3954, "B-", -1452.0, 249.074850428, 1.269);
    static const Species Es_249("Es-249", "Es", 51, 150, 99, 249, 7474.0, "B-", -2344.0, 249.076409, 32.0);
    static const Species Fm_249("Fm-249", "Fm", 49, 149, 100, 249, 7461.8649, "B-", -3661.8091, 249.078926042, 6.668);
    static const Species Md_249("Md-249", "Md", 47, 148, 101, 249, 7444.0169, "B-", -4606.0, 249.082857155, 176.516);
    static const Species No_249("No-249", "No", 45, 147, 102, 249, 7422.0, "B-", std::numeric_limits<double>::quiet_NaN(), 249.087802, 300.0);
    static const Species Cm_250("Cm-250", "Cm", 58, 154, 96, 250, 7478.9385, "B-", 37.582, 250.078357541, 11.029);
    static const Species Bk_250("Bk-250", "Bk", 56, 153, 97, 250, 7475.9594, "B-", 1781.6696, 250.078317195, 3.11);
    static const Species Cf_250("Cf-250", "Cf", 54, 152, 98, 250, 7479.9567, "B-", -2055.0, 250.076404494, 1.65);
    static const Species Es_250("Es-250", "Es", 52, 151, 99, 250, 7469.0, "B-", -847.0, 250.078611, 107.0);
    static const Species Fm_250("Fm-250", "Fm", 50, 150, 100, 250, 7462.0905, "B-", -4326.9476, 250.079519765, 8.468);
    static const Species Md_250("Md-250", "Md", 48, 149, 101, 250, 7441.6533, "B-", -3167.0, 250.084164934, 97.606);
    static const Species No_250("No-250", "No", 46, 148, 102, 250, 7426.0, "B-", std::numeric_limits<double>::quiet_NaN(), 250.087565, 215.0);
    static const Species Cm_251("Cm-251", "Cm", 59, 155, 96, 251, 7466.7233, "B-", 1420.0, 251.082284988, 24.367);
    static const Species Bk_251("Bk-251", "Bk", 57, 154, 97, 251, 7469.2637, "B-", 1093.0, 251.080760555, 11.523);
    static const Species Cf_251("Cf-251", "Cf", 55, 153, 98, 251, 7470.5014, "B-", -376.566, 251.079587171, 4.187);
    static const Species Es_251("Es-251", "Es", 53, 152, 99, 251, 7465.8842, "B-", -1447.261, 251.079991431, 5.676);
    static const Species Fm_251("Fm-251", "Fm", 51, 151, 100, 251, 7457.0013, "B-", -3007.9406, 251.08154513, 15.342);
    static const Species Md_251("Md-251", "Md", 49, 150, 101, 251, 7441.9005, "B-", -3882.0, 251.084774287, 20.31);
    static const Species No_251("No-251", "No", 47, 149, 102, 251, 7423.0, "B-", -4981.0, 251.088942, 194.0);
    static const Species Lr_251("Lr-251", "Lr", 45, 148, 103, 251, 7400.0, "B-", std::numeric_limits<double>::quiet_NaN(), 251.094289, 215.0);
    static const Species Cm_252("Cm-252", "Cm", 60, 156, 96, 252, 7460.0, "B-", 521.0, 252.08487, 320.0);
    static const Species Bk_252("Bk-252", "Bk", 58, 155, 97, 252, 7459.0, "B-", 2500.0, 252.08431, 215.0);
    static const Species Cf_252("Cf-252", "Cf", 56, 154, 98, 252, 7465.3474, "B-", -1260.0, 252.081626507, 2.531);
    static const Species Es_252("Es-252", "Es", 54, 153, 99, 252, 7457.2428, "B-", 477.9998, 252.082979173, 53.736);
    static const Species Fm_252("Fm-252", "Fm", 52, 152, 100, 252, 7456.0351, "B-", -3650.5075, 252.082466019, 5.604);
    static const Species Md_252("Md-252", "Md", 50, 151, 101, 252, 7438.4444, "B-", -2404.2523, 252.086385, 98.0);
    static const Species No_252("No-252", "No", 48, 150, 102, 252, 7425.7992, "B-", -5666.0, 252.08896607, 9.975);
    static const Species Lr_252("Lr-252", "Lr", 46, 149, 103, 252, 7400.0, "B-", std::numeric_limits<double>::quiet_NaN(), 252.095048, 198.0);
    static const Species Bk_253("Bk-253", "Bk", 59, 156, 97, 253, 7451.0, "B-", 1627.0, 253.08688, 385.0);
    static const Species Cf_253("Cf-253", "Cf", 57, 155, 98, 253, 7454.8297, "B-", 291.0753, 253.085133723, 4.57);
    static const Species Es_253("Es-253", "Es", 55, 154, 99, 253, 7452.8879, "B-", -335.0623, 253.084821241, 1.341);
    static const Species Fm_253("Fm-253", "Fm", 53, 153, 100, 253, 7448.4712, "B-", -1827.0, 253.085180945, 1.662);
    static const Species Md_253("Md-253", "Md", 51, 152, 101, 253, 7438.0, "B-", -3186.0, 253.087143, 34.0);
    static const Species No_253("No-253", "No", 49, 151, 102, 253, 7422.4719, "B-", -4164.7752, 253.09056278, 7.42);
    static const Species Lr_253("Lr-253", "Lr", 47, 150, 103, 253, 7402.918, "B-", -5118.0, 253.09503385, 176.634);
    static const Species Rf_253("Rf-253", "Rf", 45, 149, 104, 253, 7380.0, "B-", std::numeric_limits<double>::quiet_NaN(), 253.100528, 440.0);
    static const Species Bk_254("Bk-254", "Bk", 60, 157, 97, 254, 7440.0, "B-", 3052.0, 254.0906, 320.0);
    static const Species Cf_254("Cf-254", "Cf", 58, 156, 98, 254, 7449.2259, "B-", -652.7561, 254.087323575, 12.304);
    static const Species Es_254("Es-254", "Es", 56, 155, 99, 254, 7443.5759, "B-", 1091.63, 254.088024337, 3.152);
    static const Species Fm_254("Fm-254", "Fm", 54, 154, 100, 254, 7444.7936, "B-", -2550.0, 254.086852424, 1.978);
    static const Species Md_254("Md-254", "Md", 52, 153, 101, 254, 7432.0, "B-", -1271.0, 254.08959, 107.0);
    static const Species No_254("No-254", "No", 50, 152, 102, 254, 7423.5909, "B-", -4922.5753, 254.090954211, 10.367);
    static const Species Lr_254("Lr-254", "Lr", 48, 151, 103, 254, 7401.1306, "B-", -3555.0, 254.096238813, 98.026);
    static const Species Rf_254("Rf-254", "Rf", 46, 150, 104, 254, 7384.0, "B-", std::numeric_limits<double>::quiet_NaN(), 254.100055, 304.0);
    static const Species Cf_255("Cf-255", "Cf", 59, 157, 98, 255, 7438.0, "B-", 720.0, 255.091046, 215.0);
    static const Species Es_255("Es-255", "Es", 57, 156, 99, 255, 7437.8216, "B-", 288.7717, 255.090273504, 11.612);
    static const Species Fm_255("Fm-255", "Fm", 55, 155, 100, 255, 7435.886, "B-", -1041.6037, 255.089963495, 4.223);
    static const Species Md_255("Md-255", "Md", 53, 154, 101, 255, 7428.7333, "B-", -1969.8648, 255.091081702, 5.976);
    static const Species No_255("No-255", "No", 51, 153, 102, 255, 7417.9403, "B-", -3135.3716, 255.093196439, 15.079);
    static const Species Lr_255("Lr-255", "Lr", 49, 152, 103, 255, 7402.5767, "B-", -4382.0, 255.096562399, 19.0);
    static const Species Rf_255("Rf-255", "Rf", 47, 151, 104, 255, 7382.0, "B-", -5265.0, 255.101267, 194.0);
    static const Species Db_255("Db-255", "Db", 45, 150, 105, 255, 7359.0, "B-", std::numeric_limits<double>::quiet_NaN(), 255.106919, 304.0);
    static const Species Cf_256("Cf-256", "Cf", 60, 158, 98, 256, 7432.0, "B-", -144.0, 256.093442, 338.0);
    static const Species Es_256("Es-256", "Es", 58, 157, 99, 256, 7428.0, "B-", 1700.0, 256.093597, 107.0);
    static const Species Fm_256("Fm-256", "Fm", 56, 156, 100, 256, 7431.7888, "B-", -1971.0, 256.091771699, 3.241);
    static const Species Md_256("Md-256", "Md", 54, 155, 101, 256, 7421.0, "B-", -367.0, 256.093888, 133.0);
    static const Species No_256("No-256", "No", 52, 154, 102, 256, 7416.5429, "B-", -3923.5573, 256.094281912, 8.103);
    static const Species Lr_256("Lr-256", "Lr", 50, 153, 103, 256, 7398.1605, "B-", -2475.3893, 256.098494024, 89.0);
    static const Species Rf_256("Rf-256", "Rf", 48, 152, 104, 256, 7385.4349, "B-", -6076.0, 256.101151464, 19.16);
    static const Species Db_256("Db-256", "Db", 46, 151, 105, 256, 7359.0, "B-", std::numeric_limits<double>::quiet_NaN(), 256.107674, 201.0);
    static const Species Es_257("Es-257", "Es", 59, 158, 99, 257, 7422.0, "B-", 813.0, 257.095979, 441.0);
    static const Species Fm_257("Fm-257", "Fm", 57, 157, 100, 257, 7422.1942, "B-", -402.3347, 257.095105419, 4.669);
    static const Species Md_257("Md-257", "Md", 55, 156, 101, 257, 7417.5845, "B-", -1254.5923, 257.095537343, 1.683);
    static const Species No_257("No-257", "No", 53, 155, 102, 257, 7409.6587, "B-", -2418.0, 257.096884203, 6.652);
    static const Species Lr_257("Lr-257", "Lr", 51, 154, 103, 257, 7397.0, "B-", -3201.0, 257.09948, 47.0);
    static const Species Rf_257("Rf-257", "Rf", 49, 153, 104, 257, 7381.7053, "B-", -4287.8969, 257.102916796, 11.612);
    static const Species Db_257("Db-257", "Db", 47, 152, 105, 257, 7361.9767, "B-", std::numeric_limits<double>::quiet_NaN(), 257.107520042, 176.741);
    static const Species Es_258("Es-258", "Es", 60, 159, 99, 258, 7412.0, "B-", 2276.0, 258.09952, 430.0);
    static const Species Fm_258("Fm-258", "Fm", 58, 158, 100, 258, 7418.0, "B-", -1264.0, 258.097077, 215.0);
    static const Species Md_258("Md-258", "Md", 56, 157, 101, 258, 7409.6615, "B-", 213.0, 258.098433634, 3.729);
    static const Species No_258("No-258", "No", 54, 156, 102, 258, 7407.0, "B-", -3304.0, 258.098205, 107.0);
    static const Species Lr_258("Lr-258", "Lr", 52, 155, 103, 258, 7392.0, "B-", -1562.0, 258.101753, 109.0);
    static const Species Rf_258("Rf-258", "Rf", 50, 154, 104, 258, 7382.5257, "B-", -5163.3651, 258.103429895, 17.288);
    static const Species Db_258("Db-258", "Db", 48, 153, 105, 258, 7359.4803, "B-", -3788.0, 258.108972995, 98.613);
    static const Species Sg_258("Sg-258", "Sg", 46, 152, 106, 258, 7342.0, "B-", std::numeric_limits<double>::quiet_NaN(), 258.11304, 443.0);
    static const Species Fm_259("Fm-259", "Fm", 59, 159, 100, 259, 7407.0, "B-", 140.0, 259.100596, 304.0);
    static const Species Md_259("Md-259", "Md", 57, 158, 101, 259, 7405.0, "B-", -515.0, 259.100445, 108.0);
    static const Species No_259("No-259", "No", 55, 157, 102, 259, 7399.9714, "B-", -1771.0, 259.100998364, 6.829);
    static const Species Lr_259("Lr-259", "Lr", 53, 156, 103, 259, 7390.0, "B-", -2516.0, 259.1029, 76.0);
    static const Species Rf_259("Rf-259", "Rf", 51, 155, 104, 259, 7377.0, "B-", -3624.0, 259.105601, 78.0);
    static const Species Db_259("Db-259", "Db", 49, 154, 105, 259, 7360.3626, "B-", -4528.0, 259.109491859, 60.854);
    static const Species Sg_259("Sg-259", "Sg", 47, 153, 106, 259, 7340.0, "B-", std::numeric_limits<double>::quiet_NaN(), 259.114353, 194.0);
    static const Species Fm_260("Fm-260", "Fm", 60, 160, 100, 260, 7402.0, "B-", -784.0, 260.102809, 467.0);
    static const Species Md_260("Md-260", "Md", 58, 159, 101, 260, 7396.0, "B-", 940.0, 260.10365, 339.0);
    static const Species No_260("No-260", "No", 56, 158, 102, 260, 7397.0, "B-", -2667.0, 260.102641, 215.0);
    static const Species Lr_260("Lr-260", "Lr", 54, 157, 103, 260, 7383.0, "B-", -871.0, 260.105504, 134.0);
    static const Species Rf_260("Rf-260", "Rf", 52, 156, 104, 260, 7377.0, "B-", -4525.0, 260.10644, 215.0);
    static const Species Db_260("Db-260", "Db", 50, 155, 105, 260, 7357.0, "B-", -2875.0, 260.111297, 100.0);
    static const Species Sg_260("Sg-260", "Sg", 48, 154, 106, 260, 7342.5632, "B-", -6576.0, 260.114383435, 22.045);
    static const Species Bh_260("Bh-260", "Bh", 46, 153, 107, 260, 7314.0, "B-", std::numeric_limits<double>::quiet_NaN(), 260.121443, 211.0);
    static const Species Md_261("Md-261", "Md", 59, 160, 101, 261, 7391.0, "B-", 123.0, 261.105828, 546.0);
    static const Species No_261("No-261", "No", 57, 159, 102, 261, 7388.0, "B-", -1102.0, 261.105696, 215.0);
    static const Species Lr_261("Lr-261", "Lr", 55, 158, 103, 261, 7381.0, "B-", -1761.0, 261.106879, 215.0);
    static const Species Rf_261("Rf-261", "Rf", 53, 157, 104, 261, 7371.3858, "B-", -2990.0, 261.108769591, 70.492);
    static const Species Db_261("Db-261", "Db", 51, 156, 105, 261, 7357.0, "B-", -3697.0, 261.111979, 118.0);
    static const Species Sg_261("Sg-261", "Sg", 49, 155, 106, 261, 7339.771, "B-", -5074.4052, 261.115948135, 19.853);
    static const Species Bh_261("Bh-261", "Bh", 47, 154, 107, 261, 7317.3313, "B-", std::numeric_limits<double>::quiet_NaN(), 261.121395733, 193.026);
    static const Species Md_262("Md-262", "Md", 60, 161, 101, 262, 7382.0, "B-", 1566.0, 262.109144, 481.0);
    static const Species No_262("No-262", "No", 58, 160, 102, 262, 7385.0, "B-", -2004.0, 262.107463, 387.0);
    static const Species Lr_262("Lr-262", "Lr", 56, 159, 103, 262, 7374.0, "B-", -287.0, 262.109615, 215.0);
    static const Species Rf_262("Rf-262", "Rf", 54, 158, 104, 262, 7370.0, "B-", -3861.0, 262.109923, 240.0);
    static const Species Db_262("Db-262", "Db", 52, 157, 105, 262, 7352.0, "B-", -2116.0, 262.114067, 154.0);
    static const Species Sg_262("Sg-262", "Sg", 50, 156, 106, 262, 7341.1736, "B-", -5883.0463, 262.116338978, 23.797);
    static const Species Bh_262("Bh-262", "Bh", 48, 155, 107, 262, 7315.7331, "B-", std::numeric_limits<double>::quiet_NaN(), 262.122654688, 99.919);
    static const Species No_263("No-263", "No", 59, 161, 102, 263, 7376.0, "B-", -540.0, 263.110714, 526.0);
    static const Species Lr_263("Lr-263", "Lr", 57, 160, 103, 263, 7371.0, "B-", -1087.0, 263.111293, 240.0);
    static const Species Rf_263("Rf-263", "Rf", 55, 159, 104, 263, 7364.0, "B-", -2353.0, 263.112461, 164.0);
    static const Species Db_263("Db-263", "Db", 53, 158, 105, 263, 7352.0, "B-", -3085.0, 263.114987, 180.0);
    static const Species Sg_263("Sg-263", "Sg", 51, 157, 106, 263, 7337.0, "B-", -4301.0, 263.118299, 101.0);
    static const Species Bh_263("Bh-263", "Bh", 49, 156, 107, 263, 7318.0, "B-", -5182.0, 263.122916, 328.0);
    static const Species Hs_263("Hs-263", "Hs", 47, 155, 108, 263, 7295.0, "B-", std::numeric_limits<double>::quiet_NaN(), 263.128479, 212.0);
    static const Species No_264("No-264", "No", 60, 162, 102, 264, 7371.0, "B-", -1364.0, 264.112734, 634.0);
    static const Species Lr_264("Lr-264", "Lr", 58, 161, 103, 264, 7363.0, "B-", 300.0, 264.114198, 468.0);
    static const Species Rf_264("Rf-264", "Rf", 56, 160, 104, 264, 7361.0, "B-", -3187.0, 264.113876, 387.0);
    static const Species Db_264("Db-264", "Db", 54, 159, 105, 264, 7346.0, "B-", -1521.0, 264.117297, 253.0);
    static const Species Sg_264("Sg-264", "Sg", 52, 158, 106, 264, 7338.0, "B-", -5175.0, 264.11893, 304.0);
    static const Species Bh_264("Bh-264", "Bh", 50, 157, 107, 264, 7315.0, "B-", -3605.0, 264.124486, 190.0);
    static const Species Hs_264("Hs-264", "Hs", 48, 156, 108, 264, 7298.3762, "B-", std::numeric_limits<double>::quiet_NaN(), 264.12835633, 31.005);
    static const Species Lr_265("Lr-265", "Lr", 59, 162, 103, 265, 7359.0, "B-", -457.0, 265.116193, 587.0);
    static const Species Rf_265("Rf-265", "Rf", 57, 161, 104, 265, 7354.0, "B-", -1692.0, 265.116683, 387.0);
    static const Species Db_265("Db-265", "Db", 55, 160, 105, 265, 7345.0, "B-", -2412.0, 265.1185, 240.0);
    static const Species Sg_265("Sg-265", "Sg", 53, 159, 106, 265, 7333.0, "B-", -3601.0, 265.121089, 149.0);
    static const Species Bh_265("Bh-265", "Bh", 51, 158, 107, 265, 7316.0, "B-", -4505.0, 265.124955, 257.0);
    static const Species Hs_265("Hs-265", "Hs", 49, 157, 108, 265, 7296.2474, "B-", -5724.0, 265.129791744, 25.719);
    static const Species Mt_265("Mt-265", "Mt", 47, 156, 109, 265, 7272.0, "B-", std::numeric_limits<double>::quiet_NaN(), 265.135937, 471.0);
    static const Species Lr_266("Lr-266", "Lr", 60, 163, 103, 266, 7349.0, "B-", 1526.0, 266.119874, 579.0);
    static const Species Rf_266("Rf-266", "Rf", 58, 162, 104, 266, 7351.0, "B-", -2604.0, 266.118236, 443.0);
    static const Species Db_266("Db-266", "Db", 56, 161, 105, 266, 7339.0, "B-", -877.0, 266.121032, 304.0);
    static const Species Sg_266("Sg-266", "Sg", 54, 160, 106, 266, 7332.0, "B-", -4487.0, 266.121973, 263.0);
    static const Species Bh_266("Bh-266", "Bh", 52, 159, 107, 266, 7313.0, "B-", -3036.0, 266.12679, 175.0);
    static const Species Hs_266("Hs-266", "Hs", 50, 158, 108, 266, 7298.2611, "B-", -6533.0066, 266.130048783, 29.099);
    static const Species Mt_266("Mt-266", "Mt", 48, 157, 109, 266, 7270.7598, "B-", std::numeric_limits<double>::quiet_NaN(), 266.137062253, 103.568);
    static const Species Rf_267("Rf-267", "Rf", 59, 163, 104, 267, 7342.0, "B-", -570.0, 267.121787, 617.0);
    static const Species Db_267("Db-267", "Db", 57, 162, 105, 267, 7337.0, "B-", -1792.0, 267.122399, 402.0);
    static const Species Sg_267("Sg-267", "Sg", 55, 161, 106, 267, 7327.0, "B-", -2958.0, 267.124323, 281.0);
    static const Species Bh_267("Bh-267", "Bh", 53, 160, 107, 267, 7313.0, "B-", -3893.0, 267.127499, 282.0);
    static const Species Hs_267("Hs-267", "Hs", 51, 159, 108, 267, 7295.0, "B-", -5133.0, 267.131678, 102.0);
    static const Species Mt_267("Mt-267", "Mt", 49, 158, 109, 267, 7273.0, "B-", -6089.0, 267.137189, 540.0);
    static const Species Ds_267("Ds-267", "Ds", 47, 157, 110, 267, 7248.0, "B-", std::numeric_limits<double>::quiet_NaN(), 267.143726, 219.0);
    static const Species Rf_268("Rf-268", "Rf", 60, 164, 104, 268, 7337.0, "B-", -1584.0, 268.123968, 711.0);
    static const Species Db_268("Db-268", "Db", 58, 163, 105, 268, 7328.0, "B-", 260.0, 268.125669, 568.0);
    static const Species Sg_268("Sg-268", "Sg", 56, 162, 106, 268, 7326.0, "B-", -3907.0, 268.125389, 504.0);
    static const Species Bh_268("Bh-268", "Bh", 54, 161, 107, 268, 7309.0, "B-", -2261.0, 268.129584, 410.0);
    static const Species Hs_268("Hs-268", "Hs", 52, 160, 108, 268, 7297.0, "B-", -6183.0, 268.132011, 322.0);
    static const Species Mt_268("Mt-268", "Mt", 50, 159, 109, 268, 7271.0, "B-", -4497.0, 268.138649, 250.0);
    static const Species Ds_268("Ds-268", "Ds", 48, 158, 110, 268, 7252.0, "B-", std::numeric_limits<double>::quiet_NaN(), 268.143477, 324.0);
    static const Species Db_269("Db-269", "Db", 59, 164, 105, 269, 7323.0, "B-", -544.0, 269.127911, 669.0);
    static const Species Sg_269("Sg-269", "Sg", 57, 163, 106, 269, 7318.0, "B-", -1785.0, 269.128495, 395.0);
    static const Species Bh_269("Bh-269", "Bh", 55, 162, 107, 269, 7309.0, "B-", -3016.0, 269.130411, 402.0);
    static const Species Hs_269("Hs-269", "Hs", 53, 161, 108, 269, 7294.0, "B-", -4807.0, 269.133649, 141.0);
    static const Species Mt_269("Mt-269", "Mt", 51, 160, 109, 269, 7274.0, "B-", -5535.0, 269.138809, 335.0);
    static const Species Ds_269("Ds-269", "Ds", 49, 159, 110, 269, 7250.1551, "B-", std::numeric_limits<double>::quiet_NaN(), 269.144750965, 33.712);
    static const Species Db_270("Db-270", "Db", 60, 165, 105, 270, 7314.0, "B-", 966.0, 270.131399, 617.0);
    static const Species Sg_270("Sg-270", "Sg", 58, 164, 106, 270, 7314.0, "B-", -2799.0, 270.130362, 492.0);
    static const Species Bh_270("Bh-270", "Bh", 56, 163, 107, 270, 7301.0, "B-", -882.0, 270.133366, 320.0);
    static const Species Hs_270("Hs-270", "Hs", 54, 162, 108, 270, 7295.0, "B-", -5597.0, 270.134313, 266.0);
    static const Species Mt_270("Mt-270", "Mt", 52, 161, 109, 270, 7271.0, "B-", -3973.0, 270.140322, 205.0);
    static const Species Ds_270("Ds-270", "Ds", 50, 160, 110, 270, 7253.7634, "B-", std::numeric_limits<double>::quiet_NaN(), 270.14458662, 42.163);
    static const Species Sg_271("Sg-271", "Sg", 59, 165, 106, 271, 7305.0, "B-", -1242.0, 271.133782, 634.0);
    static const Species Bh_271("Bh-271", "Bh", 57, 164, 107, 271, 7298.0, "B-", -1832.0, 271.135115, 412.0);
    static const Species Hs_271("Hs-271", "Hs", 55, 163, 108, 271, 7288.0, "B-", -3409.0, 271.137082, 296.0);
    static const Species Mt_271("Mt-271", "Mt", 53, 162, 109, 271, 7273.0, "B-", -4853.0, 271.140741, 354.0);
    static const Species Ds_271("Ds-271", "Ds", 51, 161, 110, 271, 7252.0, "B-", std::numeric_limits<double>::quiet_NaN(), 271.145951, 104.0);
    static const Species Sg_272("Sg-272", "Sg", 60, 166, 106, 272, 7301.0, "B-", -2267.0, 272.135825, 743.0);
    static const Species Bh_272("Bh-272", "Bh", 58, 165, 107, 272, 7290.0, "B-", -217.0, 272.138259, 571.0);
    static const Species Hs_272("Hs-272", "Hs", 56, 164, 108, 272, 7286.0, "B-", -4477.0, 272.138492, 547.0);
    static const Species Mt_272("Mt-272", "Mt", 54, 163, 109, 272, 7267.0, "B-", -2601.0, 272.143298, 521.0);
    static const Species Ds_272("Ds-272", "Ds", 52, 162, 110, 272, 7255.0, "B-", -6690.0, 272.146091, 456.0);
    static const Species Rg_272("Rg-272", "Rg", 50, 161, 111, 272, 7227.0, "B-", std::numeric_limits<double>::quiet_NaN(), 272.153273, 251.0);
    static const Species Sg_273("Sg-273", "Sg", 61, 167, 106, 273, 7292.0, "B-", -763.0, 273.139475, 429.0);
    static const Species Bh_273("Bh-273", "Bh", 59, 166, 107, 273, 7286.0, "B-", -1084.0, 273.140294, 703.0);
    static const Species Hs_273("Hs-273", "Hs", 57, 165, 108, 273, 7279.0, "B-", -3015.0, 273.141458, 401.0);
    static const Species Mt_273("Mt-273", "Mt", 55, 164, 109, 273, 7265.0, "B-", -3503.0, 273.144695, 455.0);
    static const Species Ds_273("Ds-273", "Ds", 53, 163, 110, 273, 7250.0, "B-", -4600.0, 273.148455, 152.0);
    static const Species Rg_273("Rg-273", "Rg", 51, 162, 111, 273, 7230.0, "B-", std::numeric_limits<double>::quiet_NaN(), 273.153393, 429.0);
    static const Species Bh_274("Bh-274", "Bh", 60, 167, 107, 274, 7278.0, "B-", 356.0, 274.143599, 620.0);
    static const Species Hs_274("Hs-274", "Hs", 58, 166, 108, 274, 7276.0, "B-", -3843.0, 274.143217, 504.0);
    static const Species Mt_274("Mt-274", "Mt", 56, 165, 109, 274, 7259.0, "B-", -1948.0, 274.147343, 404.0);
    static const Species Ds_274("Ds-274", "Ds", 54, 164, 110, 274, 7249.0, "B-", -5415.0, 274.149434, 418.0);
    static const Species Rg_274("Rg-274", "Rg", 52, 163, 111, 274, 7227.0, "B-", std::numeric_limits<double>::quiet_NaN(), 274.155247, 225.0);
    static const Species Bh_275("Bh-275", "Bh", 61, 168, 107, 275, 7273.0, "B-", -712.0, 275.145766, 644.0);
    static const Species Hs_275("Hs-275", "Hs", 59, 167, 108, 275, 7268.0, "B-", -2275.0, 275.14653, 637.0);
    static const Species Mt_275("Mt-275", "Mt", 57, 166, 109, 275, 7257.0, "B-", -2899.0, 275.148972, 416.0);
    static const Species Ds_275("Ds-275", "Ds", 55, 165, 110, 275, 7243.0, "B-", -3729.0, 275.152085, 366.0);
    static const Species Rg_275("Rg-275", "Rg", 53, 164, 111, 275, 7227.0, "B-", std::numeric_limits<double>::quiet_NaN(), 275.156088, 479.0);
    static const Species Bh_276("Bh-276", "Bh", 62, 169, 107, 276, 7265.0, "B-", 765.0, 276.149169, 644.0);
    static const Species Hs_276("Hs-276", "Hs", 60, 168, 108, 276, 7265.0, "B-", -3127.0, 276.148348, 773.0);
    static const Species Mt_276("Mt-276", "Mt", 58, 167, 109, 276, 7250.0, "B-", -1227.0, 276.151705, 571.0);
    static const Species Ds_276("Ds-276", "Ds", 56, 166, 110, 276, 7243.0, "B-", -4847.0, 276.153022, 588.0);
    static const Species Rg_276("Rg-276", "Rg", 54, 165, 111, 276, 7223.0, "B-", -2974.0, 276.158226, 675.0);
    static const Species Cn_276("Cn-276", "Cn", 52, 164, 112, 276, 7209.0, "B-", std::numeric_limits<double>::quiet_NaN(), 276.161418, 537.0);
    static const Species Bh_277("Bh-277", "Bh", 63, 170, 107, 277, 7260.0, "B-", -275.0, 277.151477, 644.0);
    static const Species Hs_277("Hs-277", "Hs", 61, 169, 108, 277, 7256.0, "B-", -1633.0, 277.151772, 480.0);
    static const Species Mt_277("Mt-277", "Mt", 59, 168, 109, 277, 7247.0, "B-", -2084.0, 277.153525, 711.0);
    static const Species Ds_277("Ds-277", "Ds", 57, 167, 110, 277, 7237.0, "B-", -3315.0, 277.155763, 421.0);
    static const Species Rg_277("Rg-277", "Rg", 55, 166, 111, 277, 7222.0, "B-", -3925.0, 277.159322, 504.0);
    static const Species Cn_277("Cn-277", "Cn", 53, 165, 112, 277, 7205.0, "B-", std::numeric_limits<double>::quiet_NaN(), 277.163535, 165.0);
    static const Species Bh_278("Bh-278", "Bh", 64, 171, 107, 278, 7251.0, "B-", 1150.0, 278.154988, 429.0);
    static const Species Hs_278("Hs-278", "Hs", 62, 170, 108, 278, 7252.0, "B-", -2547.0, 278.153753, 322.0);
    static const Species Mt_278("Mt-278", "Mt", 60, 169, 109, 278, 7240.0, "B-", -484.0, 278.156487, 621.0);
    static const Species Ds_278("Ds-278", "Ds", 58, 168, 110, 278, 7236.0, "B-", -4270.0, 278.157007, 548.0);
    static const Species Rg_278("Rg-278", "Rg", 56, 167, 111, 278, 7218.0, "B-", -2321.0, 278.16159, 417.0);
    static const Species Cn_278("Cn-278", "Cn", 54, 166, 112, 278, 7206.0, "B-", -6188.0, 278.164083, 470.0);
    static const Species Nh_278("Nh-278", "Nh", 52, 165, 113, 278, 7181.0, "B-", std::numeric_limits<double>::quiet_NaN(), 278.170725, 240.0);
    static const Species Hs_279("Hs-279", "Hs", 63, 171, 108, 279, 7243.0, "B-", -1085.0, 279.157274, 644.0);
    static const Species Mt_279("Mt-279", "Mt", 61, 170, 109, 279, 7237.0, "B-", -1439.0, 279.158439, 720.0);
    static const Species Ds_279("Ds-279", "Ds", 59, 169, 110, 279, 7229.0, "B-", -2697.0, 279.159984, 649.0);
    static const Species Rg_279("Rg-279", "Rg", 57, 168, 111, 279, 7216.0, "B-", -3299.0, 279.16288, 453.0);
    static const Species Cn_279("Cn-279", "Cn", 55, 167, 112, 279, 7202.0, "B-", -4439.0, 279.166422, 424.0);
    static const Species Nh_279("Nh-279", "Nh", 53, 166, 113, 279, 7183.0, "B-", std::numeric_limits<double>::quiet_NaN(), 279.171187, 644.0);
    static const Species Hs_280("Hs-280", "Hs", 64, 172, 108, 280, 7239.0, "B-", -2090.0, 280.159335, 644.0);
    static const Species Mt_280("Mt-280", "Mt", 62, 171, 109, 280, 7229.0, "B-", 190.0, 280.161579, 644.0);
    static const Species Ds_280("Ds-280", "Ds", 60, 170, 110, 280, 7227.0, "B-", -3566.0, 280.161375, 803.0);
    static const Species Rg_280("Rg-280", "Rg", 58, 169, 111, 280, 7212.0, "B-", -1768.0, 280.165204, 571.0);
    static const Species Cn_280("Cn-280", "Cn", 56, 168, 112, 280, 7202.0, "B-", -5585.0, 280.167102, 626.0);
    static const Species Nh_280("Nh-280", "Nh", 54, 167, 113, 280, 7180.0, "B-", std::numeric_limits<double>::quiet_NaN(), 280.173098, 429.0);
    static const Species Mt_281("Mt-281", "Mt", 63, 172, 109, 281, 7225.0, "B-", -873.0, 281.163608, 644.0);
    static const Species Ds_281("Ds-281", "Ds", 61, 171, 110, 281, 7220.0, "B-", -2060.0, 281.164545, 529.0);
    static const Species Rg_281("Rg-281", "Rg", 59, 170, 111, 281, 7209.0, "B-", -2614.0, 281.166757, 831.0);
    static const Species Cn_281("Cn-281", "Cn", 57, 169, 112, 281, 7197.0, "B-", -3863.0, 281.169563, 427.0);
    static const Species Nh_281("Nh-281", "Nh", 55, 168, 113, 281, 7181.0, "B-", std::numeric_limits<double>::quiet_NaN(), 281.17371, 322.0);
    static const Species Mt_282("Mt-282", "Mt", 64, 173, 109, 282, 7218.0, "B-", 665.0, 282.166888, 480.0);
    static const Species Ds_282("Ds-282", "Ds", 62, 172, 110, 282, 7217.0, "B-", -2952.0, 282.166174, 322.0);
    static const Species Rg_282("Rg-282", "Rg", 60, 171, 111, 282, 7204.0, "B-", -1084.0, 282.169343, 631.0);
    static const Species Cn_282("Cn-282", "Cn", 58, 170, 112, 282, 7197.0, "B-", -4903.0, 282.170507, 588.0);
    static const Species Nh_282("Nh-282", "Nh", 56, 169, 113, 282, 7177.0, "B-", std::numeric_limits<double>::quiet_NaN(), 282.17577, 430.0);
    static const Species Ds_283("Ds-283", "Ds", 63, 173, 110, 283, 7210.0, "B-", -1550.0, 283.169437, 537.0);
    static const Species Rg_283("Rg-283", "Rg", 61, 172, 111, 283, 7201.0, "B-", -1957.0, 283.171101, 728.0);
    static const Species Cn_283("Cn-283", "Cn", 59, 171, 112, 283, 7192.0, "B-", -3226.0, 283.173202, 660.0);
    static const Species Nh_283("Nh-283", "Nh", 57, 170, 113, 283, 7177.0, "B-", std::numeric_limits<double>::quiet_NaN(), 283.176666, 469.0);
    static const Species Ds_284("Ds-284", "Ds", 64, 174, 110, 284, 7207.0, "B-", -2510.0, 284.171187, 537.0);
    static const Species Rg_284("Rg-284", "Rg", 62, 173, 111, 284, 7195.0, "B-", -445.0, 284.173882, 537.0);
    static const Species Cn_284("Cn-284", "Cn", 60, 172, 112, 284, 7191.0, "B-", -4176.0, 284.17436, 819.0);
    static const Species Nh_284("Nh-284", "Nh", 58, 171, 113, 284, 7173.0, "B-", -2188.0, 284.178843, 573.0);
    static const Species Fl_284("Fl-284", "Fl", 56, 170, 114, 284, 7163.0, "B-", std::numeric_limits<double>::quiet_NaN(), 284.181192, 704.0);
    static const Species Rg_285("Rg-285", "Rg", 63, 174, 111, 285, 7192.0, "B-", -1357.0, 285.175771, 644.0);
    static const Species Cn_285("Cn-285", "Cn", 61, 173, 112, 285, 7185.0, "B-", -2682.0, 285.177227, 544.0);
    static const Species Nh_285("Nh-285", "Nh", 59, 172, 113, 285, 7172.0, "B-", -3164.0, 285.180106, 832.0);
    static const Species Fl_285("Fl-285", "Fl", 57, 171, 114, 285, 7159.0, "B-", std::numeric_limits<double>::quiet_NaN(), 285.183503, 433.0);
    static const Species Rg_286("Rg-286", "Rg", 64, 175, 111, 286, 7185.0, "B-", 61.0, 286.178756, 492.0);
    static const Species Cn_286("Cn-286", "Cn", 62, 174, 112, 286, 7183.0, "B-", -3507.0, 286.178691, 751.0);
    static const Species Nh_286("Nh-286", "Nh", 60, 173, 113, 286, 7168.0, "B-", -1649.0, 286.182456, 634.0);
    static const Species Fl_286("Fl-286", "Fl", 58, 172, 114, 286, 7159.0, "B-", std::numeric_limits<double>::quiet_NaN(), 286.184226, 590.0);
    static const Species Cn_287("Cn-287", "Cn", 63, 175, 112, 287, 7176.0, "B-", -2085.0, 287.181826, 751.0);
    static const Species Nh_287("Nh-287", "Nh", 61, 174, 113, 287, 7166.0, "B-", -2474.0, 287.184064, 759.0);
    static const Species Fl_287("Fl-287", "Fl", 59, 173, 114, 287, 7155.0, "B-", -3819.0, 287.18672, 663.0);
    static const Species Mc_287("Mc-287", "Mc", 57, 172, 115, 287, 7139.0, "B-", std::numeric_limits<double>::quiet_NaN(), 287.19082, 475.0);
    static const Species Cn_288("Cn-288", "Cn", 64, 176, 112, 288, 7174.0, "B-", -3039.0, 288.183501, 751.0);
    static const Species Nh_288("Nh-288", "Nh", 62, 175, 113, 288, 7160.0, "B-", -947.0, 288.186764, 751.0);
    static const Species Fl_288("Fl-288", "Fl", 60, 174, 114, 288, 7154.0, "B-", -4749.0, 288.187781, 819.0);
    static const Species Mc_288("Mc-288", "Mc", 58, 173, 115, 288, 7135.0, "B-", std::numeric_limits<double>::quiet_NaN(), 288.192879, 575.0);
    static const Species Nh_289("Nh-289", "Nh", 63, 176, 113, 289, 7158.0, "B-", -1915.0, 289.188461, 537.0);
    static const Species Fl_289("Fl-289", "Fl", 61, 175, 114, 289, 7149.0, "B-", -3217.0, 289.190517, 548.0);
    static const Species Mc_289("Mc-289", "Mc", 59, 174, 115, 289, 7135.0, "B-", -3774.0, 289.193971, 834.0);
    static const Species Lv_289("Lv-289", "Lv", 57, 173, 116, 289, 7119.0, "B-", std::numeric_limits<double>::quiet_NaN(), 289.198023, 540.0);
    static const Species Nh_290("Nh-290", "Nh", 64, 177, 113, 290, 7152.0, "B-", -416.0, 290.191429, 503.0);
    static const Species Fl_290("Fl-290", "Fl", 62, 176, 114, 290, 7147.0, "B-", -4061.0, 290.191875, 752.0);
    static const Species Mc_290("Mc-290", "Mc", 60, 175, 115, 290, 7131.0, "B-", -2236.0, 290.196235, 635.0);
    static const Species Lv_290("Lv-290", "Lv", 58, 174, 116, 290, 7120.0, "B-", std::numeric_limits<double>::quiet_NaN(), 290.198635, 593.0);
    static const Species Fl_291("Fl-291", "Fl", 63, 177, 114, 291, 7141.0, "B-", -2680.0, 291.194848, 751.0);
    static const Species Mc_291("Mc-291", "Mc", 61, 176, 115, 291, 7129.0, "B-", -3064.0, 291.197725, 789.0);
    static const Species Lv_291("Lv-291", "Lv", 59, 175, 116, 291, 7116.0, "B-", -4409.0, 291.201014, 669.0);
    static const Species Ts_291("Ts-291", "Ts", 57, 174, 117, 291, 7098.0, "B-", std::numeric_limits<double>::quiet_NaN(), 291.205748, 640.0);
    static const Species Mc_292("Mc-292", "Mc", 62, 177, 115, 292, 7124.0, "B-", -1533.0, 292.200323, 751.0);
    static const Species Lv_292("Lv-292", "Lv", 60, 176, 116, 292, 7116.0, "B-", -5488.0, 292.201969, 819.0);
    static const Species Ts_292("Ts-292", "Ts", 58, 175, 117, 292, 7095.0, "B-", std::numeric_limits<double>::quiet_NaN(), 292.207861, 718.0);
    static const Species Lv_293("Lv-293", "Lv", 61, 177, 116, 293, 7111.0, "B-", -3860.0, 293.204583, 553.0);
    static const Species Ts_293("Ts-293", "Ts", 59, 176, 117, 293, 7095.0, "B-", -4374.0, 293.208727, 835.0);
    static const Species Og_293("Og-293", "Og", 57, 175, 118, 293, 7078.0, "B-", std::numeric_limits<double>::quiet_NaN(), 293.213423, 761.0);
    static const Species Ts_294("Ts-294", "Ts", 60, 177, 117, 294, 7092.0, "B-", -2923.0, 294.21084, 637.0);
    static const Species Og_294("Og-294", "Og", 58, 176, 118, 294, 7079.0, "B-", std::numeric_limits<double>::quiet_NaN(), 294.213979, 594.0);
    static const Species Og_295("Og-295", "Og", 59, 177, 118, 295, 7076.0, "B-", std::numeric_limits<double>::quiet_NaN(), 295.216178, 703.0);
    static const std::unordered_map<std::string, Species> species = {
        {"n-1", n_1},
        {"H-1", H_1},
        {"H-2", H_2},
        {"H-3", H_3},
        {"He-3", He_3},
        {"Li-3", Li_3},
        {"H-4", H_4},
        {"He-4", He_4},
        {"Li-4", Li_4},
        {"H-5", H_5},
        {"He-5", He_5},
        {"Li-5", Li_5},
        {"Be-5", Be_5},
        {"H-6", H_6},
        {"He-6", He_6},
        {"Li-6", Li_6},
        {"Be-6", Be_6},
        {"B-6", B_6},
        {"H-7", H_7},
        {"He-7", He_7},
        {"Li-7", Li_7},
        {"Be-7", Be_7},
        {"B-7", B_7},
        {"He-8", He_8},
        {"Li-8", Li_8},
        {"Be-8", Be_8},
        {"B-8", B_8},
        {"C-8", C_8},
        {"He-9", He_9},
        {"Li-9", Li_9},
        {"Be-9", Be_9},
        {"B-9", B_9},
        {"C-9", C_9},
        {"He-10", He_10},
        {"Li-10", Li_10},
        {"Be-10", Be_10},
        {"B-10", B_10},
        {"C-10", C_10},
        {"N-10", N_10},
        {"Li-11", Li_11},
        {"Be-11", Be_11},
        {"B-11", B_11},
        {"C-11", C_11},
        {"N-11", N_11},
        {"O-11", O_11},
        {"Li-12", Li_12},
        {"Be-12", Be_12},
        {"B-12", B_12},
        {"C-12", C_12},
        {"N-12", N_12},
        {"O-12", O_12},
        {"Li-13", Li_13},
        {"Be-13", Be_13},
        {"B-13", B_13},
        {"C-13", C_13},
        {"N-13", N_13},
        {"O-13", O_13},
        {"F-13", F_13},
        {"Be-14", Be_14},
        {"B-14", B_14},
        {"C-14", C_14},
        {"N-14", N_14},
        {"O-14", O_14},
        {"F-14", F_14},
        {"Be-15", Be_15},
        {"B-15", B_15},
        {"C-15", C_15},
        {"N-15", N_15},
        {"O-15", O_15},
        {"F-15", F_15},
        {"Ne-15", Ne_15},
        {"Be-16", Be_16},
        {"B-16", B_16},
        {"C-16", C_16},
        {"N-16", N_16},
        {"O-16", O_16},
        {"F-16", F_16},
        {"Ne-16", Ne_16},
        {"B-17", B_17},
        {"C-17", C_17},
        {"N-17", N_17},
        {"O-17", O_17},
        {"F-17", F_17},
        {"Ne-17", Ne_17},
        {"Na-17", Na_17},
        {"B-18", B_18},
        {"C-18", C_18},
        {"N-18", N_18},
        {"O-18", O_18},
        {"F-18", F_18},
        {"Ne-18", Ne_18},
        {"Na-18", Na_18},
        {"B-19", B_19},
        {"C-19", C_19},
        {"N-19", N_19},
        {"O-19", O_19},
        {"F-19", F_19},
        {"Ne-19", Ne_19},
        {"Na-19", Na_19},
        {"Mg-19", Mg_19},
        {"B-20", B_20},
        {"C-20", C_20},
        {"N-20", N_20},
        {"O-20", O_20},
        {"F-20", F_20},
        {"Ne-20", Ne_20},
        {"Na-20", Na_20},
        {"Mg-20", Mg_20},
        {"B-21", B_21},
        {"C-21", C_21},
        {"N-21", N_21},
        {"O-21", O_21},
        {"F-21", F_21},
        {"Ne-21", Ne_21},
        {"Na-21", Na_21},
        {"Mg-21", Mg_21},
        {"Al-21", Al_21},
        {"C-22", C_22},
        {"N-22", N_22},
        {"O-22", O_22},
        {"F-22", F_22},
        {"Ne-22", Ne_22},
        {"Na-22", Na_22},
        {"Mg-22", Mg_22},
        {"Al-22", Al_22},
        {"Si-22", Si_22},
        {"C-23", C_23},
        {"N-23", N_23},
        {"O-23", O_23},
        {"F-23", F_23},
        {"Ne-23", Ne_23},
        {"Na-23", Na_23},
        {"Mg-23", Mg_23},
        {"Al-23", Al_23},
        {"Si-23", Si_23},
        {"N-24", N_24},
        {"O-24", O_24},
        {"F-24", F_24},
        {"Ne-24", Ne_24},
        {"Na-24", Na_24},
        {"Mg-24", Mg_24},
        {"Al-24", Al_24},
        {"Si-24", Si_24},
        {"P-24", P_24},
        {"N-25", N_25},
        {"O-25", O_25},
        {"F-25", F_25},
        {"Ne-25", Ne_25},
        {"Na-25", Na_25},
        {"Mg-25", Mg_25},
        {"Al-25", Al_25},
        {"Si-25", Si_25},
        {"P-25", P_25},
        {"O-26", O_26},
        {"F-26", F_26},
        {"Ne-26", Ne_26},
        {"Na-26", Na_26},
        {"Mg-26", Mg_26},
        {"Al-26", Al_26},
        {"Si-26", Si_26},
        {"P-26", P_26},
        {"S-26", S_26},
        {"O-27", O_27},
        {"F-27", F_27},
        {"Ne-27", Ne_27},
        {"Na-27", Na_27},
        {"Mg-27", Mg_27},
        {"Al-27", Al_27},
        {"Si-27", Si_27},
        {"P-27", P_27},
        {"S-27", S_27},
        {"O-28", O_28},
        {"F-28", F_28},
        {"Ne-28", Ne_28},
        {"Na-28", Na_28},
        {"Mg-28", Mg_28},
        {"Al-28", Al_28},
        {"Si-28", Si_28},
        {"P-28", P_28},
        {"S-28", S_28},
        {"Cl-28", Cl_28},
        {"F-29", F_29},
        {"Ne-29", Ne_29},
        {"Na-29", Na_29},
        {"Mg-29", Mg_29},
        {"Al-29", Al_29},
        {"Si-29", Si_29},
        {"P-29", P_29},
        {"S-29", S_29},
        {"Cl-29", Cl_29},
        {"Ar-29", Ar_29},
        {"F-30", F_30},
        {"Ne-30", Ne_30},
        {"Na-30", Na_30},
        {"Mg-30", Mg_30},
        {"Al-30", Al_30},
        {"Si-30", Si_30},
        {"P-30", P_30},
        {"S-30", S_30},
        {"Cl-30", Cl_30},
        {"Ar-30", Ar_30},
        {"F-31", F_31},
        {"Ne-31", Ne_31},
        {"Na-31", Na_31},
        {"Mg-31", Mg_31},
        {"Al-31", Al_31},
        {"Si-31", Si_31},
        {"P-31", P_31},
        {"S-31", S_31},
        {"Cl-31", Cl_31},
        {"Ar-31", Ar_31},
        {"K-31", K_31},
        {"Ne-32", Ne_32},
        {"Na-32", Na_32},
        {"Mg-32", Mg_32},
        {"Al-32", Al_32},
        {"Si-32", Si_32},
        {"P-32", P_32},
        {"S-32", S_32},
        {"Cl-32", Cl_32},
        {"Ar-32", Ar_32},
        {"K-32", K_32},
        {"Ne-33", Ne_33},
        {"Na-33", Na_33},
        {"Mg-33", Mg_33},
        {"Al-33", Al_33},
        {"Si-33", Si_33},
        {"P-33", P_33},
        {"S-33", S_33},
        {"Cl-33", Cl_33},
        {"Ar-33", Ar_33},
        {"K-33", K_33},
        {"Ca-33", Ca_33},
        {"Ne-34", Ne_34},
        {"Na-34", Na_34},
        {"Mg-34", Mg_34},
        {"Al-34", Al_34},
        {"Si-34", Si_34},
        {"P-34", P_34},
        {"S-34", S_34},
        {"Cl-34", Cl_34},
        {"Ar-34", Ar_34},
        {"K-34", K_34},
        {"Ca-34", Ca_34},
        {"Na-35", Na_35},
        {"Mg-35", Mg_35},
        {"Al-35", Al_35},
        {"Si-35", Si_35},
        {"P-35", P_35},
        {"S-35", S_35},
        {"Cl-35", Cl_35},
        {"Ar-35", Ar_35},
        {"K-35", K_35},
        {"Ca-35", Ca_35},
        {"Sc-35", Sc_35},
        {"Na-36", Na_36},
        {"Mg-36", Mg_36},
        {"Al-36", Al_36},
        {"Si-36", Si_36},
        {"P-36", P_36},
        {"S-36", S_36},
        {"Cl-36", Cl_36},
        {"Ar-36", Ar_36},
        {"K-36", K_36},
        {"Ca-36", Ca_36},
        {"Sc-36", Sc_36},
        {"Na-37", Na_37},
        {"Mg-37", Mg_37},
        {"Al-37", Al_37},
        {"Si-37", Si_37},
        {"P-37", P_37},
        {"S-37", S_37},
        {"Cl-37", Cl_37},
        {"Ar-37", Ar_37},
        {"K-37", K_37},
        {"Ca-37", Ca_37},
        {"Sc-37", Sc_37},
        {"Ti-37", Ti_37},
        {"Na-38", Na_38},
        {"Mg-38", Mg_38},
        {"Al-38", Al_38},
        {"Si-38", Si_38},
        {"P-38", P_38},
        {"S-38", S_38},
        {"Cl-38", Cl_38},
        {"Ar-38", Ar_38},
        {"K-38", K_38},
        {"Ca-38", Ca_38},
        {"Sc-38", Sc_38},
        {"Ti-38", Ti_38},
        {"Na-39", Na_39},
        {"Mg-39", Mg_39},
        {"Al-39", Al_39},
        {"Si-39", Si_39},
        {"P-39", P_39},
        {"S-39", S_39},
        {"Cl-39", Cl_39},
        {"Ar-39", Ar_39},
        {"K-39", K_39},
        {"Ca-39", Ca_39},
        {"Sc-39", Sc_39},
        {"Ti-39", Ti_39},
        {"V-39", V_39},
        {"Mg-40", Mg_40},
        {"Al-40", Al_40},
        {"Si-40", Si_40},
        {"P-40", P_40},
        {"S-40", S_40},
        {"Cl-40", Cl_40},
        {"Ar-40", Ar_40},
        {"K-40", K_40},
        {"Ca-40", Ca_40},
        {"Sc-40", Sc_40},
        {"Ti-40", Ti_40},
        {"V-40", V_40},
        {"Mg-41", Mg_41},
        {"Al-41", Al_41},
        {"Si-41", Si_41},
        {"P-41", P_41},
        {"S-41", S_41},
        {"Cl-41", Cl_41},
        {"Ar-41", Ar_41},
        {"K-41", K_41},
        {"Ca-41", Ca_41},
        {"Sc-41", Sc_41},
        {"Ti-41", Ti_41},
        {"V-41", V_41},
        {"Cr-41", Cr_41},
        {"Al-42", Al_42},
        {"Si-42", Si_42},
        {"P-42", P_42},
        {"S-42", S_42},
        {"Cl-42", Cl_42},
        {"Ar-42", Ar_42},
        {"K-42", K_42},
        {"Ca-42", Ca_42},
        {"Sc-42", Sc_42},
        {"Ti-42", Ti_42},
        {"V-42", V_42},
        {"Cr-42", Cr_42},
        {"Al-43", Al_43},
        {"Si-43", Si_43},
        {"P-43", P_43},
        {"S-43", S_43},
        {"Cl-43", Cl_43},
        {"Ar-43", Ar_43},
        {"K-43", K_43},
        {"Ca-43", Ca_43},
        {"Sc-43", Sc_43},
        {"Ti-43", Ti_43},
        {"V-43", V_43},
        {"Cr-43", Cr_43},
        {"Mn-43", Mn_43},
        {"Si-44", Si_44},
        {"P-44", P_44},
        {"S-44", S_44},
        {"Cl-44", Cl_44},
        {"Ar-44", Ar_44},
        {"K-44", K_44},
        {"Ca-44", Ca_44},
        {"Sc-44", Sc_44},
        {"Ti-44", Ti_44},
        {"V-44", V_44},
        {"Cr-44", Cr_44},
        {"Mn-44", Mn_44},
        {"Si-45", Si_45},
        {"P-45", P_45},
        {"S-45", S_45},
        {"Cl-45", Cl_45},
        {"Ar-45", Ar_45},
        {"K-45", K_45},
        {"Ca-45", Ca_45},
        {"Sc-45", Sc_45},
        {"Ti-45", Ti_45},
        {"V-45", V_45},
        {"Cr-45", Cr_45},
        {"Mn-45", Mn_45},
        {"Fe-45", Fe_45},
        {"P-46", P_46},
        {"S-46", S_46},
        {"Cl-46", Cl_46},
        {"Ar-46", Ar_46},
        {"K-46", K_46},
        {"Ca-46", Ca_46},
        {"Sc-46", Sc_46},
        {"Ti-46", Ti_46},
        {"V-46", V_46},
        {"Cr-46", Cr_46},
        {"Mn-46", Mn_46},
        {"Fe-46", Fe_46},
        {"P-47", P_47},
        {"S-47", S_47},
        {"Cl-47", Cl_47},
        {"Ar-47", Ar_47},
        {"K-47", K_47},
        {"Ca-47", Ca_47},
        {"Sc-47", Sc_47},
        {"Ti-47", Ti_47},
        {"V-47", V_47},
        {"Cr-47", Cr_47},
        {"Mn-47", Mn_47},
        {"Fe-47", Fe_47},
        {"Co-47", Co_47},
        {"S-48", S_48},
        {"Cl-48", Cl_48},
        {"Ar-48", Ar_48},
        {"K-48", K_48},
        {"Ca-48", Ca_48},
        {"Sc-48", Sc_48},
        {"Ti-48", Ti_48},
        {"V-48", V_48},
        {"Cr-48", Cr_48},
        {"Mn-48", Mn_48},
        {"Fe-48", Fe_48},
        {"Co-48", Co_48},
        {"Ni-48", Ni_48},
        {"S-49", S_49},
        {"Cl-49", Cl_49},
        {"Ar-49", Ar_49},
        {"K-49", K_49},
        {"Ca-49", Ca_49},
        {"Sc-49", Sc_49},
        {"Ti-49", Ti_49},
        {"V-49", V_49},
        {"Cr-49", Cr_49},
        {"Mn-49", Mn_49},
        {"Fe-49", Fe_49},
        {"Co-49", Co_49},
        {"Ni-49", Ni_49},
        {"Cl-50", Cl_50},
        {"Ar-50", Ar_50},
        {"K-50", K_50},
        {"Ca-50", Ca_50},
        {"Sc-50", Sc_50},
        {"Ti-50", Ti_50},
        {"V-50", V_50},
        {"Cr-50", Cr_50},
        {"Mn-50", Mn_50},
        {"Fe-50", Fe_50},
        {"Co-50", Co_50},
        {"Ni-50", Ni_50},
        {"Cl-51", Cl_51},
        {"Ar-51", Ar_51},
        {"K-51", K_51},
        {"Ca-51", Ca_51},
        {"Sc-51", Sc_51},
        {"Ti-51", Ti_51},
        {"V-51", V_51},
        {"Cr-51", Cr_51},
        {"Mn-51", Mn_51},
        {"Fe-51", Fe_51},
        {"Co-51", Co_51},
        {"Ni-51", Ni_51},
        {"Cl-52", Cl_52},
        {"Ar-52", Ar_52},
        {"K-52", K_52},
        {"Ca-52", Ca_52},
        {"Sc-52", Sc_52},
        {"Ti-52", Ti_52},
        {"V-52", V_52},
        {"Cr-52", Cr_52},
        {"Mn-52", Mn_52},
        {"Fe-52", Fe_52},
        {"Co-52", Co_52},
        {"Ni-52", Ni_52},
        {"Cu-52", Cu_52},
        {"Ar-53", Ar_53},
        {"K-53", K_53},
        {"Ca-53", Ca_53},
        {"Sc-53", Sc_53},
        {"Ti-53", Ti_53},
        {"V-53", V_53},
        {"Cr-53", Cr_53},
        {"Mn-53", Mn_53},
        {"Fe-53", Fe_53},
        {"Co-53", Co_53},
        {"Ni-53", Ni_53},
        {"Cu-53", Cu_53},
        {"Ar-54", Ar_54},
        {"K-54", K_54},
        {"Ca-54", Ca_54},
        {"Sc-54", Sc_54},
        {"Ti-54", Ti_54},
        {"V-54", V_54},
        {"Cr-54", Cr_54},
        {"Mn-54", Mn_54},
        {"Fe-54", Fe_54},
        {"Co-54", Co_54},
        {"Ni-54", Ni_54},
        {"Cu-54", Cu_54},
        {"Zn-54", Zn_54},
        {"K-55", K_55},
        {"Ca-55", Ca_55},
        {"Sc-55", Sc_55},
        {"Ti-55", Ti_55},
        {"V-55", V_55},
        {"Cr-55", Cr_55},
        {"Mn-55", Mn_55},
        {"Fe-55", Fe_55},
        {"Co-55", Co_55},
        {"Ni-55", Ni_55},
        {"Cu-55", Cu_55},
        {"Zn-55", Zn_55},
        {"K-56", K_56},
        {"Ca-56", Ca_56},
        {"Sc-56", Sc_56},
        {"Ti-56", Ti_56},
        {"V-56", V_56},
        {"Cr-56", Cr_56},
        {"Mn-56", Mn_56},
        {"Fe-56", Fe_56},
        {"Co-56", Co_56},
        {"Ni-56", Ni_56},
        {"Cu-56", Cu_56},
        {"Zn-56", Zn_56},
        {"Ga-56", Ga_56},
        {"K-57", K_57},
        {"Ca-57", Ca_57},
        {"Sc-57", Sc_57},
        {"Ti-57", Ti_57},
        {"V-57", V_57},
        {"Cr-57", Cr_57},
        {"Mn-57", Mn_57},
        {"Fe-57", Fe_57},
        {"Co-57", Co_57},
        {"Ni-57", Ni_57},
        {"Cu-57", Cu_57},
        {"Zn-57", Zn_57},
        {"Ga-57", Ga_57},
        {"K-58", K_58},
        {"Ca-58", Ca_58},
        {"Sc-58", Sc_58},
        {"Ti-58", Ti_58},
        {"V-58", V_58},
        {"Cr-58", Cr_58},
        {"Mn-58", Mn_58},
        {"Fe-58", Fe_58},
        {"Co-58", Co_58},
        {"Ni-58", Ni_58},
        {"Cu-58", Cu_58},
        {"Zn-58", Zn_58},
        {"Ga-58", Ga_58},
        {"Ge-58", Ge_58},
        {"K-59", K_59},
        {"Ca-59", Ca_59},
        {"Sc-59", Sc_59},
        {"Ti-59", Ti_59},
        {"V-59", V_59},
        {"Cr-59", Cr_59},
        {"Mn-59", Mn_59},
        {"Fe-59", Fe_59},
        {"Co-59", Co_59},
        {"Ni-59", Ni_59},
        {"Cu-59", Cu_59},
        {"Zn-59", Zn_59},
        {"Ga-59", Ga_59},
        {"Ge-59", Ge_59},
        {"Ca-60", Ca_60},
        {"Sc-60", Sc_60},
        {"Ti-60", Ti_60},
        {"V-60", V_60},
        {"Cr-60", Cr_60},
        {"Mn-60", Mn_60},
        {"Fe-60", Fe_60},
        {"Co-60", Co_60},
        {"Ni-60", Ni_60},
        {"Cu-60", Cu_60},
        {"Zn-60", Zn_60},
        {"Ga-60", Ga_60},
        {"Ge-60", Ge_60},
        {"As-60", As_60},
        {"Ca-61", Ca_61},
        {"Sc-61", Sc_61},
        {"Ti-61", Ti_61},
        {"V-61", V_61},
        {"Cr-61", Cr_61},
        {"Mn-61", Mn_61},
        {"Fe-61", Fe_61},
        {"Co-61", Co_61},
        {"Ni-61", Ni_61},
        {"Cu-61", Cu_61},
        {"Zn-61", Zn_61},
        {"Ga-61", Ga_61},
        {"Ge-61", Ge_61},
        {"As-61", As_61},
        {"Sc-62", Sc_62},
        {"Ti-62", Ti_62},
        {"V-62", V_62},
        {"Cr-62", Cr_62},
        {"Mn-62", Mn_62},
        {"Fe-62", Fe_62},
        {"Co-62", Co_62},
        {"Ni-62", Ni_62},
        {"Cu-62", Cu_62},
        {"Zn-62", Zn_62},
        {"Ga-62", Ga_62},
        {"Ge-62", Ge_62},
        {"As-62", As_62},
        {"Sc-63", Sc_63},
        {"Ti-63", Ti_63},
        {"V-63", V_63},
        {"Cr-63", Cr_63},
        {"Mn-63", Mn_63},
        {"Fe-63", Fe_63},
        {"Co-63", Co_63},
        {"Ni-63", Ni_63},
        {"Cu-63", Cu_63},
        {"Zn-63", Zn_63},
        {"Ga-63", Ga_63},
        {"Ge-63", Ge_63},
        {"As-63", As_63},
        {"Se-63", Se_63},
        {"Ti-64", Ti_64},
        {"V-64", V_64},
        {"Cr-64", Cr_64},
        {"Mn-64", Mn_64},
        {"Fe-64", Fe_64},
        {"Co-64", Co_64},
        {"Ni-64", Ni_64},
        {"Cu-64", Cu_64},
        {"Zn-64", Zn_64},
        {"Ga-64", Ga_64},
        {"Ge-64", Ge_64},
        {"As-64", As_64},
        {"Se-64", Se_64},
        {"Ti-65", Ti_65},
        {"V-65", V_65},
        {"Cr-65", Cr_65},
        {"Mn-65", Mn_65},
        {"Fe-65", Fe_65},
        {"Co-65", Co_65},
        {"Ni-65", Ni_65},
        {"Cu-65", Cu_65},
        {"Zn-65", Zn_65},
        {"Ga-65", Ga_65},
        {"Ge-65", Ge_65},
        {"As-65", As_65},
        {"Se-65", Se_65},
        {"Br-65", Br_65},
        {"V-66", V_66},
        {"Cr-66", Cr_66},
        {"Mn-66", Mn_66},
        {"Fe-66", Fe_66},
        {"Co-66", Co_66},
        {"Ni-66", Ni_66},
        {"Cu-66", Cu_66},
        {"Zn-66", Zn_66},
        {"Ga-66", Ga_66},
        {"Ge-66", Ge_66},
        {"As-66", As_66},
        {"Se-66", Se_66},
        {"Br-66", Br_66},
        {"V-67", V_67},
        {"Cr-67", Cr_67},
        {"Mn-67", Mn_67},
        {"Fe-67", Fe_67},
        {"Co-67", Co_67},
        {"Ni-67", Ni_67},
        {"Cu-67", Cu_67},
        {"Zn-67", Zn_67},
        {"Ga-67", Ga_67},
        {"Ge-67", Ge_67},
        {"As-67", As_67},
        {"Se-67", Se_67},
        {"Br-67", Br_67},
        {"Kr-67", Kr_67},
        {"Cr-68", Cr_68},
        {"Mn-68", Mn_68},
        {"Fe-68", Fe_68},
        {"Co-68", Co_68},
        {"Ni-68", Ni_68},
        {"Cu-68", Cu_68},
        {"Zn-68", Zn_68},
        {"Ga-68", Ga_68},
        {"Ge-68", Ge_68},
        {"As-68", As_68},
        {"Se-68", Se_68},
        {"Br-68", Br_68},
        {"Kr-68", Kr_68},
        {"Cr-69", Cr_69},
        {"Mn-69", Mn_69},
        {"Fe-69", Fe_69},
        {"Co-69", Co_69},
        {"Ni-69", Ni_69},
        {"Cu-69", Cu_69},
        {"Zn-69", Zn_69},
        {"Ga-69", Ga_69},
        {"Ge-69", Ge_69},
        {"As-69", As_69},
        {"Se-69", Se_69},
        {"Br-69", Br_69},
        {"Kr-69", Kr_69},
        {"Cr-70", Cr_70},
        {"Mn-70", Mn_70},
        {"Fe-70", Fe_70},
        {"Co-70", Co_70},
        {"Ni-70", Ni_70},
        {"Cu-70", Cu_70},
        {"Zn-70", Zn_70},
        {"Ga-70", Ga_70},
        {"Ge-70", Ge_70},
        {"As-70", As_70},
        {"Se-70", Se_70},
        {"Br-70", Br_70},
        {"Kr-70", Kr_70},
        {"Mn-71", Mn_71},
        {"Fe-71", Fe_71},
        {"Co-71", Co_71},
        {"Ni-71", Ni_71},
        {"Cu-71", Cu_71},
        {"Zn-71", Zn_71},
        {"Ga-71", Ga_71},
        {"Ge-71", Ge_71},
        {"As-71", As_71},
        {"Se-71", Se_71},
        {"Br-71", Br_71},
        {"Kr-71", Kr_71},
        {"Rb-71", Rb_71},
        {"Mn-72", Mn_72},
        {"Fe-72", Fe_72},
        {"Co-72", Co_72},
        {"Ni-72", Ni_72},
        {"Cu-72", Cu_72},
        {"Zn-72", Zn_72},
        {"Ga-72", Ga_72},
        {"Ge-72", Ge_72},
        {"As-72", As_72},
        {"Se-72", Se_72},
        {"Br-72", Br_72},
        {"Kr-72", Kr_72},
        {"Rb-72", Rb_72},
        {"Mn-73", Mn_73},
        {"Fe-73", Fe_73},
        {"Co-73", Co_73},
        {"Ni-73", Ni_73},
        {"Cu-73", Cu_73},
        {"Zn-73", Zn_73},
        {"Ga-73", Ga_73},
        {"Ge-73", Ge_73},
        {"As-73", As_73},
        {"Se-73", Se_73},
        {"Br-73", Br_73},
        {"Kr-73", Kr_73},
        {"Rb-73", Rb_73},
        {"Sr-73", Sr_73},
        {"Fe-74", Fe_74},
        {"Co-74", Co_74},
        {"Ni-74", Ni_74},
        {"Cu-74", Cu_74},
        {"Zn-74", Zn_74},
        {"Ga-74", Ga_74},
        {"Ge-74", Ge_74},
        {"As-74", As_74},
        {"Se-74", Se_74},
        {"Br-74", Br_74},
        {"Kr-74", Kr_74},
        {"Rb-74", Rb_74},
        {"Sr-74", Sr_74},
        {"Fe-75", Fe_75},
        {"Co-75", Co_75},
        {"Ni-75", Ni_75},
        {"Cu-75", Cu_75},
        {"Zn-75", Zn_75},
        {"Ga-75", Ga_75},
        {"Ge-75", Ge_75},
        {"As-75", As_75},
        {"Se-75", Se_75},
        {"Br-75", Br_75},
        {"Kr-75", Kr_75},
        {"Rb-75", Rb_75},
        {"Sr-75", Sr_75},
        {"Y-75", Y_75},
        {"Fe-76", Fe_76},
        {"Co-76", Co_76},
        {"Ni-76", Ni_76},
        {"Cu-76", Cu_76},
        {"Zn-76", Zn_76},
        {"Ga-76", Ga_76},
        {"Ge-76", Ge_76},
        {"As-76", As_76},
        {"Se-76", Se_76},
        {"Br-76", Br_76},
        {"Kr-76", Kr_76},
        {"Rb-76", Rb_76},
        {"Sr-76", Sr_76},
        {"Y-76", Y_76},
        {"Co-77", Co_77},
        {"Ni-77", Ni_77},
        {"Cu-77", Cu_77},
        {"Zn-77", Zn_77},
        {"Ga-77", Ga_77},
        {"Ge-77", Ge_77},
        {"As-77", As_77},
        {"Se-77", Se_77},
        {"Br-77", Br_77},
        {"Kr-77", Kr_77},
        {"Rb-77", Rb_77},
        {"Sr-77", Sr_77},
        {"Y-77", Y_77},
        {"Zr-77", Zr_77},
        {"Co-78", Co_78},
        {"Ni-78", Ni_78},
        {"Cu-78", Cu_78},
        {"Zn-78", Zn_78},
        {"Ga-78", Ga_78},
        {"Ge-78", Ge_78},
        {"As-78", As_78},
        {"Se-78", Se_78},
        {"Br-78", Br_78},
        {"Kr-78", Kr_78},
        {"Rb-78", Rb_78},
        {"Sr-78", Sr_78},
        {"Y-78", Y_78},
        {"Zr-78", Zr_78},
        {"Ni-79", Ni_79},
        {"Cu-79", Cu_79},
        {"Zn-79", Zn_79},
        {"Ga-79", Ga_79},
        {"Ge-79", Ge_79},
        {"As-79", As_79},
        {"Se-79", Se_79},
        {"Br-79", Br_79},
        {"Kr-79", Kr_79},
        {"Rb-79", Rb_79},
        {"Sr-79", Sr_79},
        {"Y-79", Y_79},
        {"Zr-79", Zr_79},
        {"Nb-79", Nb_79},
        {"Ni-80", Ni_80},
        {"Cu-80", Cu_80},
        {"Zn-80", Zn_80},
        {"Ga-80", Ga_80},
        {"Ge-80", Ge_80},
        {"As-80", As_80},
        {"Se-80", Se_80},
        {"Br-80", Br_80},
        {"Kr-80", Kr_80},
        {"Rb-80", Rb_80},
        {"Sr-80", Sr_80},
        {"Y-80", Y_80},
        {"Zr-80", Zr_80},
        {"Nb-80", Nb_80},
        {"Ni-81", Ni_81},
        {"Cu-81", Cu_81},
        {"Zn-81", Zn_81},
        {"Ga-81", Ga_81},
        {"Ge-81", Ge_81},
        {"As-81", As_81},
        {"Se-81", Se_81},
        {"Br-81", Br_81},
        {"Kr-81", Kr_81},
        {"Rb-81", Rb_81},
        {"Sr-81", Sr_81},
        {"Y-81", Y_81},
        {"Zr-81", Zr_81},
        {"Nb-81", Nb_81},
        {"Mo-81", Mo_81},
        {"Ni-82", Ni_82},
        {"Cu-82", Cu_82},
        {"Zn-82", Zn_82},
        {"Ga-82", Ga_82},
        {"Ge-82", Ge_82},
        {"As-82", As_82},
        {"Se-82", Se_82},
        {"Br-82", Br_82},
        {"Kr-82", Kr_82},
        {"Rb-82", Rb_82},
        {"Sr-82", Sr_82},
        {"Y-82", Y_82},
        {"Zr-82", Zr_82},
        {"Nb-82", Nb_82},
        {"Mo-82", Mo_82},
        {"Cu-83", Cu_83},
        {"Zn-83", Zn_83},
        {"Ga-83", Ga_83},
        {"Ge-83", Ge_83},
        {"As-83", As_83},
        {"Se-83", Se_83},
        {"Br-83", Br_83},
        {"Kr-83", Kr_83},
        {"Rb-83", Rb_83},
        {"Sr-83", Sr_83},
        {"Y-83", Y_83},
        {"Zr-83", Zr_83},
        {"Nb-83", Nb_83},
        {"Mo-83", Mo_83},
        {"Tc-83", Tc_83},
        {"Cu-84", Cu_84},
        {"Zn-84", Zn_84},
        {"Ga-84", Ga_84},
        {"Ge-84", Ge_84},
        {"As-84", As_84},
        {"Se-84", Se_84},
        {"Br-84", Br_84},
        {"Kr-84", Kr_84},
        {"Rb-84", Rb_84},
        {"Sr-84", Sr_84},
        {"Y-84", Y_84},
        {"Zr-84", Zr_84},
        {"Nb-84", Nb_84},
        {"Mo-84", Mo_84},
        {"Tc-84", Tc_84},
        {"Zn-85", Zn_85},
        {"Ga-85", Ga_85},
        {"Ge-85", Ge_85},
        {"As-85", As_85},
        {"Se-85", Se_85},
        {"Br-85", Br_85},
        {"Kr-85", Kr_85},
        {"Rb-85", Rb_85},
        {"Sr-85", Sr_85},
        {"Y-85", Y_85},
        {"Zr-85", Zr_85},
        {"Nb-85", Nb_85},
        {"Mo-85", Mo_85},
        {"Tc-85", Tc_85},
        {"Ru-85", Ru_85},
        {"Zn-86", Zn_86},
        {"Ga-86", Ga_86},
        {"Ge-86", Ge_86},
        {"As-86", As_86},
        {"Se-86", Se_86},
        {"Br-86", Br_86},
        {"Kr-86", Kr_86},
        {"Rb-86", Rb_86},
        {"Sr-86", Sr_86},
        {"Y-86", Y_86},
        {"Zr-86", Zr_86},
        {"Nb-86", Nb_86},
        {"Mo-86", Mo_86},
        {"Tc-86", Tc_86},
        {"Ru-86", Ru_86},
        {"Ga-87", Ga_87},
        {"Ge-87", Ge_87},
        {"As-87", As_87},
        {"Se-87", Se_87},
        {"Br-87", Br_87},
        {"Kr-87", Kr_87},
        {"Rb-87", Rb_87},
        {"Sr-87", Sr_87},
        {"Y-87", Y_87},
        {"Zr-87", Zr_87},
        {"Nb-87", Nb_87},
        {"Mo-87", Mo_87},
        {"Tc-87", Tc_87},
        {"Ru-87", Ru_87},
        {"Ga-88", Ga_88},
        {"Ge-88", Ge_88},
        {"As-88", As_88},
        {"Se-88", Se_88},
        {"Br-88", Br_88},
        {"Kr-88", Kr_88},
        {"Rb-88", Rb_88},
        {"Sr-88", Sr_88},
        {"Y-88", Y_88},
        {"Zr-88", Zr_88},
        {"Nb-88", Nb_88},
        {"Mo-88", Mo_88},
        {"Tc-88", Tc_88},
        {"Ru-88", Ru_88},
        {"Rh-88", Rh_88},
        {"Ge-89", Ge_89},
        {"As-89", As_89},
        {"Se-89", Se_89},
        {"Br-89", Br_89},
        {"Kr-89", Kr_89},
        {"Rb-89", Rb_89},
        {"Sr-89", Sr_89},
        {"Y-89", Y_89},
        {"Zr-89", Zr_89},
        {"Nb-89", Nb_89},
        {"Mo-89", Mo_89},
        {"Tc-89", Tc_89},
        {"Ru-89", Ru_89},
        {"Rh-89", Rh_89},
        {"Ge-90", Ge_90},
        {"As-90", As_90},
        {"Se-90", Se_90},
        {"Br-90", Br_90},
        {"Kr-90", Kr_90},
        {"Rb-90", Rb_90},
        {"Sr-90", Sr_90},
        {"Y-90", Y_90},
        {"Zr-90", Zr_90},
        {"Nb-90", Nb_90},
        {"Mo-90", Mo_90},
        {"Tc-90", Tc_90},
        {"Ru-90", Ru_90},
        {"Rh-90", Rh_90},
        {"Pd-90", Pd_90},
        {"As-91", As_91},
        {"Se-91", Se_91},
        {"Br-91", Br_91},
        {"Kr-91", Kr_91},
        {"Rb-91", Rb_91},
        {"Sr-91", Sr_91},
        {"Y-91", Y_91},
        {"Zr-91", Zr_91},
        {"Nb-91", Nb_91},
        {"Mo-91", Mo_91},
        {"Tc-91", Tc_91},
        {"Ru-91", Ru_91},
        {"Rh-91", Rh_91},
        {"Pd-91", Pd_91},
        {"As-92", As_92},
        {"Se-92", Se_92},
        {"Br-92", Br_92},
        {"Kr-92", Kr_92},
        {"Rb-92", Rb_92},
        {"Sr-92", Sr_92},
        {"Y-92", Y_92},
        {"Zr-92", Zr_92},
        {"Nb-92", Nb_92},
        {"Mo-92", Mo_92},
        {"Tc-92", Tc_92},
        {"Ru-92", Ru_92},
        {"Rh-92", Rh_92},
        {"Pd-92", Pd_92},
        {"Ag-92", Ag_92},
        {"Se-93", Se_93},
        {"Br-93", Br_93},
        {"Kr-93", Kr_93},
        {"Rb-93", Rb_93},
        {"Sr-93", Sr_93},
        {"Y-93", Y_93},
        {"Zr-93", Zr_93},
        {"Nb-93", Nb_93},
        {"Mo-93", Mo_93},
        {"Tc-93", Tc_93},
        {"Ru-93", Ru_93},
        {"Rh-93", Rh_93},
        {"Pd-93", Pd_93},
        {"Ag-93", Ag_93},
        {"Se-94", Se_94},
        {"Br-94", Br_94},
        {"Kr-94", Kr_94},
        {"Rb-94", Rb_94},
        {"Sr-94", Sr_94},
        {"Y-94", Y_94},
        {"Zr-94", Zr_94},
        {"Nb-94", Nb_94},
        {"Mo-94", Mo_94},
        {"Tc-94", Tc_94},
        {"Ru-94", Ru_94},
        {"Rh-94", Rh_94},
        {"Pd-94", Pd_94},
        {"Ag-94", Ag_94},
        {"Cd-94", Cd_94},
        {"Se-95", Se_95},
        {"Br-95", Br_95},
        {"Kr-95", Kr_95},
        {"Rb-95", Rb_95},
        {"Sr-95", Sr_95},
        {"Y-95", Y_95},
        {"Zr-95", Zr_95},
        {"Nb-95", Nb_95},
        {"Mo-95", Mo_95},
        {"Tc-95", Tc_95},
        {"Ru-95", Ru_95},
        {"Rh-95", Rh_95},
        {"Pd-95", Pd_95},
        {"Ag-95", Ag_95},
        {"Cd-95", Cd_95},
        {"Br-96", Br_96},
        {"Kr-96", Kr_96},
        {"Rb-96", Rb_96},
        {"Sr-96", Sr_96},
        {"Y-96", Y_96},
        {"Zr-96", Zr_96},
        {"Nb-96", Nb_96},
        {"Mo-96", Mo_96},
        {"Tc-96", Tc_96},
        {"Ru-96", Ru_96},
        {"Rh-96", Rh_96},
        {"Pd-96", Pd_96},
        {"Ag-96", Ag_96},
        {"Cd-96", Cd_96},
        {"In-96", In_96},
        {"Br-97", Br_97},
        {"Kr-97", Kr_97},
        {"Rb-97", Rb_97},
        {"Sr-97", Sr_97},
        {"Y-97", Y_97},
        {"Zr-97", Zr_97},
        {"Nb-97", Nb_97},
        {"Mo-97", Mo_97},
        {"Tc-97", Tc_97},
        {"Ru-97", Ru_97},
        {"Rh-97", Rh_97},
        {"Pd-97", Pd_97},
        {"Ag-97", Ag_97},
        {"Cd-97", Cd_97},
        {"In-97", In_97},
        {"Br-98", Br_98},
        {"Kr-98", Kr_98},
        {"Rb-98", Rb_98},
        {"Sr-98", Sr_98},
        {"Y-98", Y_98},
        {"Zr-98", Zr_98},
        {"Nb-98", Nb_98},
        {"Mo-98", Mo_98},
        {"Tc-98", Tc_98},
        {"Ru-98", Ru_98},
        {"Rh-98", Rh_98},
        {"Pd-98", Pd_98},
        {"Ag-98", Ag_98},
        {"Cd-98", Cd_98},
        {"In-98", In_98},
        {"Kr-99", Kr_99},
        {"Rb-99", Rb_99},
        {"Sr-99", Sr_99},
        {"Y-99", Y_99},
        {"Zr-99", Zr_99},
        {"Nb-99", Nb_99},
        {"Mo-99", Mo_99},
        {"Tc-99", Tc_99},
        {"Ru-99", Ru_99},
        {"Rh-99", Rh_99},
        {"Pd-99", Pd_99},
        {"Ag-99", Ag_99},
        {"Cd-99", Cd_99},
        {"In-99", In_99},
        {"Sn-99", Sn_99},
        {"Kr-100", Kr_100},
        {"Rb-100", Rb_100},
        {"Sr-100", Sr_100},
        {"Y-100", Y_100},
        {"Zr-100", Zr_100},
        {"Nb-100", Nb_100},
        {"Mo-100", Mo_100},
        {"Tc-100", Tc_100},
        {"Ru-100", Ru_100},
        {"Rh-100", Rh_100},
        {"Pd-100", Pd_100},
        {"Ag-100", Ag_100},
        {"Cd-100", Cd_100},
        {"In-100", In_100},
        {"Sn-100", Sn_100},
        {"Kr-101", Kr_101},
        {"Rb-101", Rb_101},
        {"Sr-101", Sr_101},
        {"Y-101", Y_101},
        {"Zr-101", Zr_101},
        {"Nb-101", Nb_101},
        {"Mo-101", Mo_101},
        {"Tc-101", Tc_101},
        {"Ru-101", Ru_101},
        {"Rh-101", Rh_101},
        {"Pd-101", Pd_101},
        {"Ag-101", Ag_101},
        {"Cd-101", Cd_101},
        {"In-101", In_101},
        {"Sn-101", Sn_101},
        {"Rb-102", Rb_102},
        {"Sr-102", Sr_102},
        {"Y-102", Y_102},
        {"Zr-102", Zr_102},
        {"Nb-102", Nb_102},
        {"Mo-102", Mo_102},
        {"Tc-102", Tc_102},
        {"Ru-102", Ru_102},
        {"Rh-102", Rh_102},
        {"Pd-102", Pd_102},
        {"Ag-102", Ag_102},
        {"Cd-102", Cd_102},
        {"In-102", In_102},
        {"Sn-102", Sn_102},
        {"Sb-102", Sb_102},
        {"Rb-103", Rb_103},
        {"Sr-103", Sr_103},
        {"Y-103", Y_103},
        {"Zr-103", Zr_103},
        {"Nb-103", Nb_103},
        {"Mo-103", Mo_103},
        {"Tc-103", Tc_103},
        {"Ru-103", Ru_103},
        {"Rh-103", Rh_103},
        {"Pd-103", Pd_103},
        {"Ag-103", Ag_103},
        {"Cd-103", Cd_103},
        {"In-103", In_103},
        {"Sn-103", Sn_103},
        {"Sb-103", Sb_103},
        {"Rb-104", Rb_104},
        {"Sr-104", Sr_104},
        {"Y-104", Y_104},
        {"Zr-104", Zr_104},
        {"Nb-104", Nb_104},
        {"Mo-104", Mo_104},
        {"Tc-104", Tc_104},
        {"Ru-104", Ru_104},
        {"Rh-104", Rh_104},
        {"Pd-104", Pd_104},
        {"Ag-104", Ag_104},
        {"Cd-104", Cd_104},
        {"In-104", In_104},
        {"Sn-104", Sn_104},
        {"Sb-104", Sb_104},
        {"Te-104", Te_104},
        {"Sr-105", Sr_105},
        {"Y-105", Y_105},
        {"Zr-105", Zr_105},
        {"Nb-105", Nb_105},
        {"Mo-105", Mo_105},
        {"Tc-105", Tc_105},
        {"Ru-105", Ru_105},
        {"Rh-105", Rh_105},
        {"Pd-105", Pd_105},
        {"Ag-105", Ag_105},
        {"Cd-105", Cd_105},
        {"In-105", In_105},
        {"Sn-105", Sn_105},
        {"Sb-105", Sb_105},
        {"Te-105", Te_105},
        {"Sr-106", Sr_106},
        {"Y-106", Y_106},
        {"Zr-106", Zr_106},
        {"Nb-106", Nb_106},
        {"Mo-106", Mo_106},
        {"Tc-106", Tc_106},
        {"Ru-106", Ru_106},
        {"Rh-106", Rh_106},
        {"Pd-106", Pd_106},
        {"Ag-106", Ag_106},
        {"Cd-106", Cd_106},
        {"In-106", In_106},
        {"Sn-106", Sn_106},
        {"Sb-106", Sb_106},
        {"Te-106", Te_106},
        {"I-106", I_106},
        {"Sr-107", Sr_107},
        {"Y-107", Y_107},
        {"Zr-107", Zr_107},
        {"Nb-107", Nb_107},
        {"Mo-107", Mo_107},
        {"Tc-107", Tc_107},
        {"Ru-107", Ru_107},
        {"Rh-107", Rh_107},
        {"Pd-107", Pd_107},
        {"Ag-107", Ag_107},
        {"Cd-107", Cd_107},
        {"In-107", In_107},
        {"Sn-107", Sn_107},
        {"Sb-107", Sb_107},
        {"Te-107", Te_107},
        {"I-107", I_107},
        {"Y-108", Y_108},
        {"Zr-108", Zr_108},
        {"Nb-108", Nb_108},
        {"Mo-108", Mo_108},
        {"Tc-108", Tc_108},
        {"Ru-108", Ru_108},
        {"Rh-108", Rh_108},
        {"Pd-108", Pd_108},
        {"Ag-108", Ag_108},
        {"Cd-108", Cd_108},
        {"In-108", In_108},
        {"Sn-108", Sn_108},
        {"Sb-108", Sb_108},
        {"Te-108", Te_108},
        {"I-108", I_108},
        {"Xe-108", Xe_108},
        {"Y-109", Y_109},
        {"Zr-109", Zr_109},
        {"Nb-109", Nb_109},
        {"Mo-109", Mo_109},
        {"Tc-109", Tc_109},
        {"Ru-109", Ru_109},
        {"Rh-109", Rh_109},
        {"Pd-109", Pd_109},
        {"Ag-109", Ag_109},
        {"Cd-109", Cd_109},
        {"In-109", In_109},
        {"Sn-109", Sn_109},
        {"Sb-109", Sb_109},
        {"Te-109", Te_109},
        {"I-109", I_109},
        {"Xe-109", Xe_109},
        {"Zr-110", Zr_110},
        {"Nb-110", Nb_110},
        {"Mo-110", Mo_110},
        {"Tc-110", Tc_110},
        {"Ru-110", Ru_110},
        {"Rh-110", Rh_110},
        {"Pd-110", Pd_110},
        {"Ag-110", Ag_110},
        {"Cd-110", Cd_110},
        {"In-110", In_110},
        {"Sn-110", Sn_110},
        {"Sb-110", Sb_110},
        {"Te-110", Te_110},
        {"I-110", I_110},
        {"Xe-110", Xe_110},
        {"Zr-111", Zr_111},
        {"Nb-111", Nb_111},
        {"Mo-111", Mo_111},
        {"Tc-111", Tc_111},
        {"Ru-111", Ru_111},
        {"Rh-111", Rh_111},
        {"Pd-111", Pd_111},
        {"Ag-111", Ag_111},
        {"Cd-111", Cd_111},
        {"In-111", In_111},
        {"Sn-111", Sn_111},
        {"Sb-111", Sb_111},
        {"Te-111", Te_111},
        {"I-111", I_111},
        {"Xe-111", Xe_111},
        {"Cs-111", Cs_111},
        {"Zr-112", Zr_112},
        {"Nb-112", Nb_112},
        {"Mo-112", Mo_112},
        {"Tc-112", Tc_112},
        {"Ru-112", Ru_112},
        {"Rh-112", Rh_112},
        {"Pd-112", Pd_112},
        {"Ag-112", Ag_112},
        {"Cd-112", Cd_112},
        {"In-112", In_112},
        {"Sn-112", Sn_112},
        {"Sb-112", Sb_112},
        {"Te-112", Te_112},
        {"I-112", I_112},
        {"Xe-112", Xe_112},
        {"Cs-112", Cs_112},
        {"Zr-113", Zr_113},
        {"Nb-113", Nb_113},
        {"Mo-113", Mo_113},
        {"Tc-113", Tc_113},
        {"Ru-113", Ru_113},
        {"Rh-113", Rh_113},
        {"Pd-113", Pd_113},
        {"Ag-113", Ag_113},
        {"Cd-113", Cd_113},
        {"In-113", In_113},
        {"Sn-113", Sn_113},
        {"Sb-113", Sb_113},
        {"Te-113", Te_113},
        {"I-113", I_113},
        {"Xe-113", Xe_113},
        {"Cs-113", Cs_113},
        {"Ba-113", Ba_113},
        {"Nb-114", Nb_114},
        {"Mo-114", Mo_114},
        {"Tc-114", Tc_114},
        {"Ru-114", Ru_114},
        {"Rh-114", Rh_114},
        {"Pd-114", Pd_114},
        {"Ag-114", Ag_114},
        {"Cd-114", Cd_114},
        {"In-114", In_114},
        {"Sn-114", Sn_114},
        {"Sb-114", Sb_114},
        {"Te-114", Te_114},
        {"I-114", I_114},
        {"Xe-114", Xe_114},
        {"Cs-114", Cs_114},
        {"Ba-114", Ba_114},
        {"Nb-115", Nb_115},
        {"Mo-115", Mo_115},
        {"Tc-115", Tc_115},
        {"Ru-115", Ru_115},
        {"Rh-115", Rh_115},
        {"Pd-115", Pd_115},
        {"Ag-115", Ag_115},
        {"Cd-115", Cd_115},
        {"In-115", In_115},
        {"Sn-115", Sn_115},
        {"Sb-115", Sb_115},
        {"Te-115", Te_115},
        {"I-115", I_115},
        {"Xe-115", Xe_115},
        {"Cs-115", Cs_115},
        {"Ba-115", Ba_115},
        {"Nb-116", Nb_116},
        {"Mo-116", Mo_116},
        {"Tc-116", Tc_116},
        {"Ru-116", Ru_116},
        {"Rh-116", Rh_116},
        {"Pd-116", Pd_116},
        {"Ag-116", Ag_116},
        {"Cd-116", Cd_116},
        {"In-116", In_116},
        {"Sn-116", Sn_116},
        {"Sb-116", Sb_116},
        {"Te-116", Te_116},
        {"I-116", I_116},
        {"Xe-116", Xe_116},
        {"Cs-116", Cs_116},
        {"Ba-116", Ba_116},
        {"La-116", La_116},
        {"Mo-117", Mo_117},
        {"Tc-117", Tc_117},
        {"Ru-117", Ru_117},
        {"Rh-117", Rh_117},
        {"Pd-117", Pd_117},
        {"Ag-117", Ag_117},
        {"Cd-117", Cd_117},
        {"In-117", In_117},
        {"Sn-117", Sn_117},
        {"Sb-117", Sb_117},
        {"Te-117", Te_117},
        {"I-117", I_117},
        {"Xe-117", Xe_117},
        {"Cs-117", Cs_117},
        {"Ba-117", Ba_117},
        {"La-117", La_117},
        {"Mo-118", Mo_118},
        {"Tc-118", Tc_118},
        {"Ru-118", Ru_118},
        {"Rh-118", Rh_118},
        {"Pd-118", Pd_118},
        {"Ag-118", Ag_118},
        {"Cd-118", Cd_118},
        {"In-118", In_118},
        {"Sn-118", Sn_118},
        {"Sb-118", Sb_118},
        {"Te-118", Te_118},
        {"I-118", I_118},
        {"Xe-118", Xe_118},
        {"Cs-118", Cs_118},
        {"Ba-118", Ba_118},
        {"La-118", La_118},
        {"Mo-119", Mo_119},
        {"Tc-119", Tc_119},
        {"Ru-119", Ru_119},
        {"Rh-119", Rh_119},
        {"Pd-119", Pd_119},
        {"Ag-119", Ag_119},
        {"Cd-119", Cd_119},
        {"In-119", In_119},
        {"Sn-119", Sn_119},
        {"Sb-119", Sb_119},
        {"Te-119", Te_119},
        {"I-119", I_119},
        {"Xe-119", Xe_119},
        {"Cs-119", Cs_119},
        {"Ba-119", Ba_119},
        {"La-119", La_119},
        {"Ce-119", Ce_119},
        {"Tc-120", Tc_120},
        {"Ru-120", Ru_120},
        {"Rh-120", Rh_120},
        {"Pd-120", Pd_120},
        {"Ag-120", Ag_120},
        {"Cd-120", Cd_120},
        {"In-120", In_120},
        {"Sn-120", Sn_120},
        {"Sb-120", Sb_120},
        {"Te-120", Te_120},
        {"I-120", I_120},
        {"Xe-120", Xe_120},
        {"Cs-120", Cs_120},
        {"Ba-120", Ba_120},
        {"La-120", La_120},
        {"Ce-120", Ce_120},
        {"Tc-121", Tc_121},
        {"Ru-121", Ru_121},
        {"Rh-121", Rh_121},
        {"Pd-121", Pd_121},
        {"Ag-121", Ag_121},
        {"Cd-121", Cd_121},
        {"In-121", In_121},
        {"Sn-121", Sn_121},
        {"Sb-121", Sb_121},
        {"Te-121", Te_121},
        {"I-121", I_121},
        {"Xe-121", Xe_121},
        {"Cs-121", Cs_121},
        {"Ba-121", Ba_121},
        {"La-121", La_121},
        {"Ce-121", Ce_121},
        {"Pr-121", Pr_121},
        {"Tc-122", Tc_122},
        {"Ru-122", Ru_122},
        {"Rh-122", Rh_122},
        {"Pd-122", Pd_122},
        {"Ag-122", Ag_122},
        {"Cd-122", Cd_122},
        {"In-122", In_122},
        {"Sn-122", Sn_122},
        {"Sb-122", Sb_122},
        {"Te-122", Te_122},
        {"I-122", I_122},
        {"Xe-122", Xe_122},
        {"Cs-122", Cs_122},
        {"Ba-122", Ba_122},
        {"La-122", La_122},
        {"Ce-122", Ce_122},
        {"Pr-122", Pr_122},
        {"Ru-123", Ru_123},
        {"Rh-123", Rh_123},
        {"Pd-123", Pd_123},
        {"Ag-123", Ag_123},
        {"Cd-123", Cd_123},
        {"In-123", In_123},
        {"Sn-123", Sn_123},
        {"Sb-123", Sb_123},
        {"Te-123", Te_123},
        {"I-123", I_123},
        {"Xe-123", Xe_123},
        {"Cs-123", Cs_123},
        {"Ba-123", Ba_123},
        {"La-123", La_123},
        {"Ce-123", Ce_123},
        {"Pr-123", Pr_123},
        {"Ru-124", Ru_124},
        {"Rh-124", Rh_124},
        {"Pd-124", Pd_124},
        {"Ag-124", Ag_124},
        {"Cd-124", Cd_124},
        {"In-124", In_124},
        {"Sn-124", Sn_124},
        {"Sb-124", Sb_124},
        {"Te-124", Te_124},
        {"I-124", I_124},
        {"Xe-124", Xe_124},
        {"Cs-124", Cs_124},
        {"Ba-124", Ba_124},
        {"La-124", La_124},
        {"Ce-124", Ce_124},
        {"Pr-124", Pr_124},
        {"Nd-124", Nd_124},
        {"Ru-125", Ru_125},
        {"Rh-125", Rh_125},
        {"Pd-125", Pd_125},
        {"Ag-125", Ag_125},
        {"Cd-125", Cd_125},
        {"In-125", In_125},
        {"Sn-125", Sn_125},
        {"Sb-125", Sb_125},
        {"Te-125", Te_125},
        {"I-125", I_125},
        {"Xe-125", Xe_125},
        {"Cs-125", Cs_125},
        {"Ba-125", Ba_125},
        {"La-125", La_125},
        {"Ce-125", Ce_125},
        {"Pr-125", Pr_125},
        {"Nd-125", Nd_125},
        {"Rh-126", Rh_126},
        {"Pd-126", Pd_126},
        {"Ag-126", Ag_126},
        {"Cd-126", Cd_126},
        {"In-126", In_126},
        {"Sn-126", Sn_126},
        {"Sb-126", Sb_126},
        {"Te-126", Te_126},
        {"I-126", I_126},
        {"Xe-126", Xe_126},
        {"Cs-126", Cs_126},
        {"Ba-126", Ba_126},
        {"La-126", La_126},
        {"Ce-126", Ce_126},
        {"Pr-126", Pr_126},
        {"Nd-126", Nd_126},
        {"Pm-126", Pm_126},
        {"Rh-127", Rh_127},
        {"Pd-127", Pd_127},
        {"Ag-127", Ag_127},
        {"Cd-127", Cd_127},
        {"In-127", In_127},
        {"Sn-127", Sn_127},
        {"Sb-127", Sb_127},
        {"Te-127", Te_127},
        {"I-127", I_127},
        {"Xe-127", Xe_127},
        {"Cs-127", Cs_127},
        {"Ba-127", Ba_127},
        {"La-127", La_127},
        {"Ce-127", Ce_127},
        {"Pr-127", Pr_127},
        {"Nd-127", Nd_127},
        {"Pm-127", Pm_127},
        {"Rh-128", Rh_128},
        {"Pd-128", Pd_128},
        {"Ag-128", Ag_128},
        {"Cd-128", Cd_128},
        {"In-128", In_128},
        {"Sn-128", Sn_128},
        {"Sb-128", Sb_128},
        {"Te-128", Te_128},
        {"I-128", I_128},
        {"Xe-128", Xe_128},
        {"Cs-128", Cs_128},
        {"Ba-128", Ba_128},
        {"La-128", La_128},
        {"Ce-128", Ce_128},
        {"Pr-128", Pr_128},
        {"Nd-128", Nd_128},
        {"Pm-128", Pm_128},
        {"Sm-128", Sm_128},
        {"Pd-129", Pd_129},
        {"Ag-129", Ag_129},
        {"Cd-129", Cd_129},
        {"In-129", In_129},
        {"Sn-129", Sn_129},
        {"Sb-129", Sb_129},
        {"Te-129", Te_129},
        {"I-129", I_129},
        {"Xe-129", Xe_129},
        {"Cs-129", Cs_129},
        {"Ba-129", Ba_129},
        {"La-129", La_129},
        {"Ce-129", Ce_129},
        {"Pr-129", Pr_129},
        {"Nd-129", Nd_129},
        {"Pm-129", Pm_129},
        {"Sm-129", Sm_129},
        {"Pd-130", Pd_130},
        {"Ag-130", Ag_130},
        {"Cd-130", Cd_130},
        {"In-130", In_130},
        {"Sn-130", Sn_130},
        {"Sb-130", Sb_130},
        {"Te-130", Te_130},
        {"I-130", I_130},
        {"Xe-130", Xe_130},
        {"Cs-130", Cs_130},
        {"Ba-130", Ba_130},
        {"La-130", La_130},
        {"Ce-130", Ce_130},
        {"Pr-130", Pr_130},
        {"Nd-130", Nd_130},
        {"Pm-130", Pm_130},
        {"Sm-130", Sm_130},
        {"Eu-130", Eu_130},
        {"Pd-131", Pd_131},
        {"Ag-131", Ag_131},
        {"Cd-131", Cd_131},
        {"In-131", In_131},
        {"Sn-131", Sn_131},
        {"Sb-131", Sb_131},
        {"Te-131", Te_131},
        {"I-131", I_131},
        {"Xe-131", Xe_131},
        {"Cs-131", Cs_131},
        {"Ba-131", Ba_131},
        {"La-131", La_131},
        {"Ce-131", Ce_131},
        {"Pr-131", Pr_131},
        {"Nd-131", Nd_131},
        {"Pm-131", Pm_131},
        {"Sm-131", Sm_131},
        {"Eu-131", Eu_131},
        {"Ag-132", Ag_132},
        {"Cd-132", Cd_132},
        {"In-132", In_132},
        {"Sn-132", Sn_132},
        {"Sb-132", Sb_132},
        {"Te-132", Te_132},
        {"I-132", I_132},
        {"Xe-132", Xe_132},
        {"Cs-132", Cs_132},
        {"Ba-132", Ba_132},
        {"La-132", La_132},
        {"Ce-132", Ce_132},
        {"Pr-132", Pr_132},
        {"Nd-132", Nd_132},
        {"Pm-132", Pm_132},
        {"Sm-132", Sm_132},
        {"Eu-132", Eu_132},
        {"Ag-133", Ag_133},
        {"Cd-133", Cd_133},
        {"In-133", In_133},
        {"Sn-133", Sn_133},
        {"Sb-133", Sb_133},
        {"Te-133", Te_133},
        {"I-133", I_133},
        {"Xe-133", Xe_133},
        {"Cs-133", Cs_133},
        {"Ba-133", Ba_133},
        {"La-133", La_133},
        {"Ce-133", Ce_133},
        {"Pr-133", Pr_133},
        {"Nd-133", Nd_133},
        {"Pm-133", Pm_133},
        {"Sm-133", Sm_133},
        {"Eu-133", Eu_133},
        {"Gd-133", Gd_133},
        {"Cd-134", Cd_134},
        {"In-134", In_134},
        {"Sn-134", Sn_134},
        {"Sb-134", Sb_134},
        {"Te-134", Te_134},
        {"I-134", I_134},
        {"Xe-134", Xe_134},
        {"Cs-134", Cs_134},
        {"Ba-134", Ba_134},
        {"La-134", La_134},
        {"Ce-134", Ce_134},
        {"Pr-134", Pr_134},
        {"Nd-134", Nd_134},
        {"Pm-134", Pm_134},
        {"Sm-134", Sm_134},
        {"Eu-134", Eu_134},
        {"Gd-134", Gd_134},
        {"Cd-135", Cd_135},
        {"In-135", In_135},
        {"Sn-135", Sn_135},
        {"Sb-135", Sb_135},
        {"Te-135", Te_135},
        {"I-135", I_135},
        {"Xe-135", Xe_135},
        {"Cs-135", Cs_135},
        {"Ba-135", Ba_135},
        {"La-135", La_135},
        {"Ce-135", Ce_135},
        {"Pr-135", Pr_135},
        {"Nd-135", Nd_135},
        {"Pm-135", Pm_135},
        {"Sm-135", Sm_135},
        {"Eu-135", Eu_135},
        {"Gd-135", Gd_135},
        {"Tb-135", Tb_135},
        {"In-136", In_136},
        {"Sn-136", Sn_136},
        {"Sb-136", Sb_136},
        {"Te-136", Te_136},
        {"I-136", I_136},
        {"Xe-136", Xe_136},
        {"Cs-136", Cs_136},
        {"Ba-136", Ba_136},
        {"La-136", La_136},
        {"Ce-136", Ce_136},
        {"Pr-136", Pr_136},
        {"Nd-136", Nd_136},
        {"Pm-136", Pm_136},
        {"Sm-136", Sm_136},
        {"Eu-136", Eu_136},
        {"Gd-136", Gd_136},
        {"Tb-136", Tb_136},
        {"In-137", In_137},
        {"Sn-137", Sn_137},
        {"Sb-137", Sb_137},
        {"Te-137", Te_137},
        {"I-137", I_137},
        {"Xe-137", Xe_137},
        {"Cs-137", Cs_137},
        {"Ba-137", Ba_137},
        {"La-137", La_137},
        {"Ce-137", Ce_137},
        {"Pr-137", Pr_137},
        {"Nd-137", Nd_137},
        {"Pm-137", Pm_137},
        {"Sm-137", Sm_137},
        {"Eu-137", Eu_137},
        {"Gd-137", Gd_137},
        {"Tb-137", Tb_137},
        {"Sn-138", Sn_138},
        {"Sb-138", Sb_138},
        {"Te-138", Te_138},
        {"I-138", I_138},
        {"Xe-138", Xe_138},
        {"Cs-138", Cs_138},
        {"Ba-138", Ba_138},
        {"La-138", La_138},
        {"Ce-138", Ce_138},
        {"Pr-138", Pr_138},
        {"Nd-138", Nd_138},
        {"Pm-138", Pm_138},
        {"Sm-138", Sm_138},
        {"Eu-138", Eu_138},
        {"Gd-138", Gd_138},
        {"Tb-138", Tb_138},
        {"Dy-138", Dy_138},
        {"Sn-139", Sn_139},
        {"Sb-139", Sb_139},
        {"Te-139", Te_139},
        {"I-139", I_139},
        {"Xe-139", Xe_139},
        {"Cs-139", Cs_139},
        {"Ba-139", Ba_139},
        {"La-139", La_139},
        {"Ce-139", Ce_139},
        {"Pr-139", Pr_139},
        {"Nd-139", Nd_139},
        {"Pm-139", Pm_139},
        {"Sm-139", Sm_139},
        {"Eu-139", Eu_139},
        {"Gd-139", Gd_139},
        {"Tb-139", Tb_139},
        {"Dy-139", Dy_139},
        {"Sn-140", Sn_140},
        {"Sb-140", Sb_140},
        {"Te-140", Te_140},
        {"I-140", I_140},
        {"Xe-140", Xe_140},
        {"Cs-140", Cs_140},
        {"Ba-140", Ba_140},
        {"La-140", La_140},
        {"Ce-140", Ce_140},
        {"Pr-140", Pr_140},
        {"Nd-140", Nd_140},
        {"Pm-140", Pm_140},
        {"Sm-140", Sm_140},
        {"Eu-140", Eu_140},
        {"Gd-140", Gd_140},
        {"Tb-140", Tb_140},
        {"Dy-140", Dy_140},
        {"Ho-140", Ho_140},
        {"Sb-141", Sb_141},
        {"Te-141", Te_141},
        {"I-141", I_141},
        {"Xe-141", Xe_141},
        {"Cs-141", Cs_141},
        {"Ba-141", Ba_141},
        {"La-141", La_141},
        {"Ce-141", Ce_141},
        {"Pr-141", Pr_141},
        {"Nd-141", Nd_141},
        {"Pm-141", Pm_141},
        {"Sm-141", Sm_141},
        {"Eu-141", Eu_141},
        {"Gd-141", Gd_141},
        {"Tb-141", Tb_141},
        {"Dy-141", Dy_141},
        {"Ho-141", Ho_141},
        {"Sb-142", Sb_142},
        {"Te-142", Te_142},
        {"I-142", I_142},
        {"Xe-142", Xe_142},
        {"Cs-142", Cs_142},
        {"Ba-142", Ba_142},
        {"La-142", La_142},
        {"Ce-142", Ce_142},
        {"Pr-142", Pr_142},
        {"Nd-142", Nd_142},
        {"Pm-142", Pm_142},
        {"Sm-142", Sm_142},
        {"Eu-142", Eu_142},
        {"Gd-142", Gd_142},
        {"Tb-142", Tb_142},
        {"Dy-142", Dy_142},
        {"Ho-142", Ho_142},
        {"Er-142", Er_142},
        {"Te-143", Te_143},
        {"I-143", I_143},
        {"Xe-143", Xe_143},
        {"Cs-143", Cs_143},
        {"Ba-143", Ba_143},
        {"La-143", La_143},
        {"Ce-143", Ce_143},
        {"Pr-143", Pr_143},
        {"Nd-143", Nd_143},
        {"Pm-143", Pm_143},
        {"Sm-143", Sm_143},
        {"Eu-143", Eu_143},
        {"Gd-143", Gd_143},
        {"Tb-143", Tb_143},
        {"Dy-143", Dy_143},
        {"Ho-143", Ho_143},
        {"Er-143", Er_143},
        {"Te-144", Te_144},
        {"I-144", I_144},
        {"Xe-144", Xe_144},
        {"Cs-144", Cs_144},
        {"Ba-144", Ba_144},
        {"La-144", La_144},
        {"Ce-144", Ce_144},
        {"Pr-144", Pr_144},
        {"Nd-144", Nd_144},
        {"Pm-144", Pm_144},
        {"Sm-144", Sm_144},
        {"Eu-144", Eu_144},
        {"Gd-144", Gd_144},
        {"Tb-144", Tb_144},
        {"Dy-144", Dy_144},
        {"Ho-144", Ho_144},
        {"Er-144", Er_144},
        {"Tm-144", Tm_144},
        {"Te-145", Te_145},
        {"I-145", I_145},
        {"Xe-145", Xe_145},
        {"Cs-145", Cs_145},
        {"Ba-145", Ba_145},
        {"La-145", La_145},
        {"Ce-145", Ce_145},
        {"Pr-145", Pr_145},
        {"Nd-145", Nd_145},
        {"Pm-145", Pm_145},
        {"Sm-145", Sm_145},
        {"Eu-145", Eu_145},
        {"Gd-145", Gd_145},
        {"Tb-145", Tb_145},
        {"Dy-145", Dy_145},
        {"Ho-145", Ho_145},
        {"Er-145", Er_145},
        {"Tm-145", Tm_145},
        {"I-146", I_146},
        {"Xe-146", Xe_146},
        {"Cs-146", Cs_146},
        {"Ba-146", Ba_146},
        {"La-146", La_146},
        {"Ce-146", Ce_146},
        {"Pr-146", Pr_146},
        {"Nd-146", Nd_146},
        {"Pm-146", Pm_146},
        {"Sm-146", Sm_146},
        {"Eu-146", Eu_146},
        {"Gd-146", Gd_146},
        {"Tb-146", Tb_146},
        {"Dy-146", Dy_146},
        {"Ho-146", Ho_146},
        {"Er-146", Er_146},
        {"Tm-146", Tm_146},
        {"I-147", I_147},
        {"Xe-147", Xe_147},
        {"Cs-147", Cs_147},
        {"Ba-147", Ba_147},
        {"La-147", La_147},
        {"Ce-147", Ce_147},
        {"Pr-147", Pr_147},
        {"Nd-147", Nd_147},
        {"Pm-147", Pm_147},
        {"Sm-147", Sm_147},
        {"Eu-147", Eu_147},
        {"Gd-147", Gd_147},
        {"Tb-147", Tb_147},
        {"Dy-147", Dy_147},
        {"Ho-147", Ho_147},
        {"Er-147", Er_147},
        {"Tm-147", Tm_147},
        {"Xe-148", Xe_148},
        {"Cs-148", Cs_148},
        {"Ba-148", Ba_148},
        {"La-148", La_148},
        {"Ce-148", Ce_148},
        {"Pr-148", Pr_148},
        {"Nd-148", Nd_148},
        {"Pm-148", Pm_148},
        {"Sm-148", Sm_148},
        {"Eu-148", Eu_148},
        {"Gd-148", Gd_148},
        {"Tb-148", Tb_148},
        {"Dy-148", Dy_148},
        {"Ho-148", Ho_148},
        {"Er-148", Er_148},
        {"Tm-148", Tm_148},
        {"Yb-148", Yb_148},
        {"Xe-149", Xe_149},
        {"Cs-149", Cs_149},
        {"Ba-149", Ba_149},
        {"La-149", La_149},
        {"Ce-149", Ce_149},
        {"Pr-149", Pr_149},
        {"Nd-149", Nd_149},
        {"Pm-149", Pm_149},
        {"Sm-149", Sm_149},
        {"Eu-149", Eu_149},
        {"Gd-149", Gd_149},
        {"Tb-149", Tb_149},
        {"Dy-149", Dy_149},
        {"Ho-149", Ho_149},
        {"Er-149", Er_149},
        {"Tm-149", Tm_149},
        {"Yb-149", Yb_149},
        {"Xe-150", Xe_150},
        {"Cs-150", Cs_150},
        {"Ba-150", Ba_150},
        {"La-150", La_150},
        {"Ce-150", Ce_150},
        {"Pr-150", Pr_150},
        {"Nd-150", Nd_150},
        {"Pm-150", Pm_150},
        {"Sm-150", Sm_150},
        {"Eu-150", Eu_150},
        {"Gd-150", Gd_150},
        {"Tb-150", Tb_150},
        {"Dy-150", Dy_150},
        {"Ho-150", Ho_150},
        {"Er-150", Er_150},
        {"Tm-150", Tm_150},
        {"Yb-150", Yb_150},
        {"Lu-150", Lu_150},
        {"Cs-151", Cs_151},
        {"Ba-151", Ba_151},
        {"La-151", La_151},
        {"Ce-151", Ce_151},
        {"Pr-151", Pr_151},
        {"Nd-151", Nd_151},
        {"Pm-151", Pm_151},
        {"Sm-151", Sm_151},
        {"Eu-151", Eu_151},
        {"Gd-151", Gd_151},
        {"Tb-151", Tb_151},
        {"Dy-151", Dy_151},
        {"Ho-151", Ho_151},
        {"Er-151", Er_151},
        {"Tm-151", Tm_151},
        {"Yb-151", Yb_151},
        {"Lu-151", Lu_151},
        {"Cs-152", Cs_152},
        {"Ba-152", Ba_152},
        {"La-152", La_152},
        {"Ce-152", Ce_152},
        {"Pr-152", Pr_152},
        {"Nd-152", Nd_152},
        {"Pm-152", Pm_152},
        {"Sm-152", Sm_152},
        {"Eu-152", Eu_152},
        {"Gd-152", Gd_152},
        {"Tb-152", Tb_152},
        {"Dy-152", Dy_152},
        {"Ho-152", Ho_152},
        {"Er-152", Er_152},
        {"Tm-152", Tm_152},
        {"Yb-152", Yb_152},
        {"Lu-152", Lu_152},
        {"Ba-153", Ba_153},
        {"La-153", La_153},
        {"Ce-153", Ce_153},
        {"Pr-153", Pr_153},
        {"Nd-153", Nd_153},
        {"Pm-153", Pm_153},
        {"Sm-153", Sm_153},
        {"Eu-153", Eu_153},
        {"Gd-153", Gd_153},
        {"Tb-153", Tb_153},
        {"Dy-153", Dy_153},
        {"Ho-153", Ho_153},
        {"Er-153", Er_153},
        {"Tm-153", Tm_153},
        {"Yb-153", Yb_153},
        {"Lu-153", Lu_153},
        {"Hf-153", Hf_153},
        {"Ba-154", Ba_154},
        {"La-154", La_154},
        {"Ce-154", Ce_154},
        {"Pr-154", Pr_154},
        {"Nd-154", Nd_154},
        {"Pm-154", Pm_154},
        {"Sm-154", Sm_154},
        {"Eu-154", Eu_154},
        {"Gd-154", Gd_154},
        {"Tb-154", Tb_154},
        {"Dy-154", Dy_154},
        {"Ho-154", Ho_154},
        {"Er-154", Er_154},
        {"Tm-154", Tm_154},
        {"Yb-154", Yb_154},
        {"Lu-154", Lu_154},
        {"Hf-154", Hf_154},
        {"La-155", La_155},
        {"Ce-155", Ce_155},
        {"Pr-155", Pr_155},
        {"Nd-155", Nd_155},
        {"Pm-155", Pm_155},
        {"Sm-155", Sm_155},
        {"Eu-155", Eu_155},
        {"Gd-155", Gd_155},
        {"Tb-155", Tb_155},
        {"Dy-155", Dy_155},
        {"Ho-155", Ho_155},
        {"Er-155", Er_155},
        {"Tm-155", Tm_155},
        {"Yb-155", Yb_155},
        {"Lu-155", Lu_155},
        {"Hf-155", Hf_155},
        {"Ta-155", Ta_155},
        {"La-156", La_156},
        {"Ce-156", Ce_156},
        {"Pr-156", Pr_156},
        {"Nd-156", Nd_156},
        {"Pm-156", Pm_156},
        {"Sm-156", Sm_156},
        {"Eu-156", Eu_156},
        {"Gd-156", Gd_156},
        {"Tb-156", Tb_156},
        {"Dy-156", Dy_156},
        {"Ho-156", Ho_156},
        {"Er-156", Er_156},
        {"Tm-156", Tm_156},
        {"Yb-156", Yb_156},
        {"Lu-156", Lu_156},
        {"Hf-156", Hf_156},
        {"Ta-156", Ta_156},
        {"La-157", La_157},
        {"Ce-157", Ce_157},
        {"Pr-157", Pr_157},
        {"Nd-157", Nd_157},
        {"Pm-157", Pm_157},
        {"Sm-157", Sm_157},
        {"Eu-157", Eu_157},
        {"Gd-157", Gd_157},
        {"Tb-157", Tb_157},
        {"Dy-157", Dy_157},
        {"Ho-157", Ho_157},
        {"Er-157", Er_157},
        {"Tm-157", Tm_157},
        {"Yb-157", Yb_157},
        {"Lu-157", Lu_157},
        {"Hf-157", Hf_157},
        {"Ta-157", Ta_157},
        {"W-157", W_157},
        {"Ce-158", Ce_158},
        {"Pr-158", Pr_158},
        {"Nd-158", Nd_158},
        {"Pm-158", Pm_158},
        {"Sm-158", Sm_158},
        {"Eu-158", Eu_158},
        {"Gd-158", Gd_158},
        {"Tb-158", Tb_158},
        {"Dy-158", Dy_158},
        {"Ho-158", Ho_158},
        {"Er-158", Er_158},
        {"Tm-158", Tm_158},
        {"Yb-158", Yb_158},
        {"Lu-158", Lu_158},
        {"Hf-158", Hf_158},
        {"Ta-158", Ta_158},
        {"W-158", W_158},
        {"Ce-159", Ce_159},
        {"Pr-159", Pr_159},
        {"Nd-159", Nd_159},
        {"Pm-159", Pm_159},
        {"Sm-159", Sm_159},
        {"Eu-159", Eu_159},
        {"Gd-159", Gd_159},
        {"Tb-159", Tb_159},
        {"Dy-159", Dy_159},
        {"Ho-159", Ho_159},
        {"Er-159", Er_159},
        {"Tm-159", Tm_159},
        {"Yb-159", Yb_159},
        {"Lu-159", Lu_159},
        {"Hf-159", Hf_159},
        {"Ta-159", Ta_159},
        {"W-159", W_159},
        {"Re-159", Re_159},
        {"Pr-160", Pr_160},
        {"Nd-160", Nd_160},
        {"Pm-160", Pm_160},
        {"Sm-160", Sm_160},
        {"Eu-160", Eu_160},
        {"Gd-160", Gd_160},
        {"Tb-160", Tb_160},
        {"Dy-160", Dy_160},
        {"Ho-160", Ho_160},
        {"Er-160", Er_160},
        {"Tm-160", Tm_160},
        {"Yb-160", Yb_160},
        {"Lu-160", Lu_160},
        {"Hf-160", Hf_160},
        {"Ta-160", Ta_160},
        {"W-160", W_160},
        {"Re-160", Re_160},
        {"Pr-161", Pr_161},
        {"Nd-161", Nd_161},
        {"Pm-161", Pm_161},
        {"Sm-161", Sm_161},
        {"Eu-161", Eu_161},
        {"Gd-161", Gd_161},
        {"Tb-161", Tb_161},
        {"Dy-161", Dy_161},
        {"Ho-161", Ho_161},
        {"Er-161", Er_161},
        {"Tm-161", Tm_161},
        {"Yb-161", Yb_161},
        {"Lu-161", Lu_161},
        {"Hf-161", Hf_161},
        {"Ta-161", Ta_161},
        {"W-161", W_161},
        {"Re-161", Re_161},
        {"Os-161", Os_161},
        {"Nd-162", Nd_162},
        {"Pm-162", Pm_162},
        {"Sm-162", Sm_162},
        {"Eu-162", Eu_162},
        {"Gd-162", Gd_162},
        {"Tb-162", Tb_162},
        {"Dy-162", Dy_162},
        {"Ho-162", Ho_162},
        {"Er-162", Er_162},
        {"Tm-162", Tm_162},
        {"Yb-162", Yb_162},
        {"Lu-162", Lu_162},
        {"Hf-162", Hf_162},
        {"Ta-162", Ta_162},
        {"W-162", W_162},
        {"Re-162", Re_162},
        {"Os-162", Os_162},
        {"Nd-163", Nd_163},
        {"Pm-163", Pm_163},
        {"Sm-163", Sm_163},
        {"Eu-163", Eu_163},
        {"Gd-163", Gd_163},
        {"Tb-163", Tb_163},
        {"Dy-163", Dy_163},
        {"Ho-163", Ho_163},
        {"Er-163", Er_163},
        {"Tm-163", Tm_163},
        {"Yb-163", Yb_163},
        {"Lu-163", Lu_163},
        {"Hf-163", Hf_163},
        {"Ta-163", Ta_163},
        {"W-163", W_163},
        {"Re-163", Re_163},
        {"Os-163", Os_163},
        {"Ir-163", Ir_163},
        {"Pm-164", Pm_164},
        {"Sm-164", Sm_164},
        {"Eu-164", Eu_164},
        {"Gd-164", Gd_164},
        {"Tb-164", Tb_164},
        {"Dy-164", Dy_164},
        {"Ho-164", Ho_164},
        {"Er-164", Er_164},
        {"Tm-164", Tm_164},
        {"Yb-164", Yb_164},
        {"Lu-164", Lu_164},
        {"Hf-164", Hf_164},
        {"Ta-164", Ta_164},
        {"W-164", W_164},
        {"Re-164", Re_164},
        {"Os-164", Os_164},
        {"Ir-164", Ir_164},
        {"Pm-165", Pm_165},
        {"Sm-165", Sm_165},
        {"Eu-165", Eu_165},
        {"Gd-165", Gd_165},
        {"Tb-165", Tb_165},
        {"Dy-165", Dy_165},
        {"Ho-165", Ho_165},
        {"Er-165", Er_165},
        {"Tm-165", Tm_165},
        {"Yb-165", Yb_165},
        {"Lu-165", Lu_165},
        {"Hf-165", Hf_165},
        {"Ta-165", Ta_165},
        {"W-165", W_165},
        {"Re-165", Re_165},
        {"Os-165", Os_165},
        {"Ir-165", Ir_165},
        {"Pt-165", Pt_165},
        {"Sm-166", Sm_166},
        {"Eu-166", Eu_166},
        {"Gd-166", Gd_166},
        {"Tb-166", Tb_166},
        {"Dy-166", Dy_166},
        {"Ho-166", Ho_166},
        {"Er-166", Er_166},
        {"Tm-166", Tm_166},
        {"Yb-166", Yb_166},
        {"Lu-166", Lu_166},
        {"Hf-166", Hf_166},
        {"Ta-166", Ta_166},
        {"W-166", W_166},
        {"Re-166", Re_166},
        {"Os-166", Os_166},
        {"Ir-166", Ir_166},
        {"Pt-166", Pt_166},
        {"Sm-167", Sm_167},
        {"Eu-167", Eu_167},
        {"Gd-167", Gd_167},
        {"Tb-167", Tb_167},
        {"Dy-167", Dy_167},
        {"Ho-167", Ho_167},
        {"Er-167", Er_167},
        {"Tm-167", Tm_167},
        {"Yb-167", Yb_167},
        {"Lu-167", Lu_167},
        {"Hf-167", Hf_167},
        {"Ta-167", Ta_167},
        {"W-167", W_167},
        {"Re-167", Re_167},
        {"Os-167", Os_167},
        {"Ir-167", Ir_167},
        {"Pt-167", Pt_167},
        {"Sm-168", Sm_168},
        {"Eu-168", Eu_168},
        {"Gd-168", Gd_168},
        {"Tb-168", Tb_168},
        {"Dy-168", Dy_168},
        {"Ho-168", Ho_168},
        {"Er-168", Er_168},
        {"Tm-168", Tm_168},
        {"Yb-168", Yb_168},
        {"Lu-168", Lu_168},
        {"Hf-168", Hf_168},
        {"Ta-168", Ta_168},
        {"W-168", W_168},
        {"Re-168", Re_168},
        {"Os-168", Os_168},
        {"Ir-168", Ir_168},
        {"Pt-168", Pt_168},
        {"Au-168", Au_168},
        {"Eu-169", Eu_169},
        {"Gd-169", Gd_169},
        {"Tb-169", Tb_169},
        {"Dy-169", Dy_169},
        {"Ho-169", Ho_169},
        {"Er-169", Er_169},
        {"Tm-169", Tm_169},
        {"Yb-169", Yb_169},
        {"Lu-169", Lu_169},
        {"Hf-169", Hf_169},
        {"Ta-169", Ta_169},
        {"W-169", W_169},
        {"Re-169", Re_169},
        {"Os-169", Os_169},
        {"Ir-169", Ir_169},
        {"Pt-169", Pt_169},
        {"Au-169", Au_169},
        {"Eu-170", Eu_170},
        {"Gd-170", Gd_170},
        {"Tb-170", Tb_170},
        {"Dy-170", Dy_170},
        {"Ho-170", Ho_170},
        {"Er-170", Er_170},
        {"Tm-170", Tm_170},
        {"Yb-170", Yb_170},
        {"Lu-170", Lu_170},
        {"Hf-170", Hf_170},
        {"Ta-170", Ta_170},
        {"W-170", W_170},
        {"Re-170", Re_170},
        {"Os-170", Os_170},
        {"Ir-170", Ir_170},
        {"Pt-170", Pt_170},
        {"Au-170", Au_170},
        {"Hg-170", Hg_170},
        {"Gd-171", Gd_171},
        {"Tb-171", Tb_171},
        {"Dy-171", Dy_171},
        {"Ho-171", Ho_171},
        {"Er-171", Er_171},
        {"Tm-171", Tm_171},
        {"Yb-171", Yb_171},
        {"Lu-171", Lu_171},
        {"Hf-171", Hf_171},
        {"Ta-171", Ta_171},
        {"W-171", W_171},
        {"Re-171", Re_171},
        {"Os-171", Os_171},
        {"Ir-171", Ir_171},
        {"Pt-171", Pt_171},
        {"Au-171", Au_171},
        {"Hg-171", Hg_171},
        {"Gd-172", Gd_172},
        {"Tb-172", Tb_172},
        {"Dy-172", Dy_172},
        {"Ho-172", Ho_172},
        {"Er-172", Er_172},
        {"Tm-172", Tm_172},
        {"Yb-172", Yb_172},
        {"Lu-172", Lu_172},
        {"Hf-172", Hf_172},
        {"Ta-172", Ta_172},
        {"W-172", W_172},
        {"Re-172", Re_172},
        {"Os-172", Os_172},
        {"Ir-172", Ir_172},
        {"Pt-172", Pt_172},
        {"Au-172", Au_172},
        {"Hg-172", Hg_172},
        {"Tb-173", Tb_173},
        {"Dy-173", Dy_173},
        {"Ho-173", Ho_173},
        {"Er-173", Er_173},
        {"Tm-173", Tm_173},
        {"Yb-173", Yb_173},
        {"Lu-173", Lu_173},
        {"Hf-173", Hf_173},
        {"Ta-173", Ta_173},
        {"W-173", W_173},
        {"Re-173", Re_173},
        {"Os-173", Os_173},
        {"Ir-173", Ir_173},
        {"Pt-173", Pt_173},
        {"Au-173", Au_173},
        {"Hg-173", Hg_173},
        {"Tb-174", Tb_174},
        {"Dy-174", Dy_174},
        {"Ho-174", Ho_174},
        {"Er-174", Er_174},
        {"Tm-174", Tm_174},
        {"Yb-174", Yb_174},
        {"Lu-174", Lu_174},
        {"Hf-174", Hf_174},
        {"Ta-174", Ta_174},
        {"W-174", W_174},
        {"Re-174", Re_174},
        {"Os-174", Os_174},
        {"Ir-174", Ir_174},
        {"Pt-174", Pt_174},
        {"Au-174", Au_174},
        {"Hg-174", Hg_174},
        {"Dy-175", Dy_175},
        {"Ho-175", Ho_175},
        {"Er-175", Er_175},
        {"Tm-175", Tm_175},
        {"Yb-175", Yb_175},
        {"Lu-175", Lu_175},
        {"Hf-175", Hf_175},
        {"Ta-175", Ta_175},
        {"W-175", W_175},
        {"Re-175", Re_175},
        {"Os-175", Os_175},
        {"Ir-175", Ir_175},
        {"Pt-175", Pt_175},
        {"Au-175", Au_175},
        {"Hg-175", Hg_175},
        {"Dy-176", Dy_176},
        {"Ho-176", Ho_176},
        {"Er-176", Er_176},
        {"Tm-176", Tm_176},
        {"Yb-176", Yb_176},
        {"Lu-176", Lu_176},
        {"Hf-176", Hf_176},
        {"Ta-176", Ta_176},
        {"W-176", W_176},
        {"Re-176", Re_176},
        {"Os-176", Os_176},
        {"Ir-176", Ir_176},
        {"Pt-176", Pt_176},
        {"Au-176", Au_176},
        {"Hg-176", Hg_176},
        {"Tl-176", Tl_176},
        {"Ho-177", Ho_177},
        {"Er-177", Er_177},
        {"Tm-177", Tm_177},
        {"Yb-177", Yb_177},
        {"Lu-177", Lu_177},
        {"Hf-177", Hf_177},
        {"Ta-177", Ta_177},
        {"W-177", W_177},
        {"Re-177", Re_177},
        {"Os-177", Os_177},
        {"Ir-177", Ir_177},
        {"Pt-177", Pt_177},
        {"Au-177", Au_177},
        {"Hg-177", Hg_177},
        {"Tl-177", Tl_177},
        {"Ho-178", Ho_178},
        {"Er-178", Er_178},
        {"Tm-178", Tm_178},
        {"Yb-178", Yb_178},
        {"Lu-178", Lu_178},
        {"Hf-178", Hf_178},
        {"Ta-178", Ta_178},
        {"W-178", W_178},
        {"Re-178", Re_178},
        {"Os-178", Os_178},
        {"Ir-178", Ir_178},
        {"Pt-178", Pt_178},
        {"Au-178", Au_178},
        {"Hg-178", Hg_178},
        {"Tl-178", Tl_178},
        {"Pb-178", Pb_178},
        {"Er-179", Er_179},
        {"Tm-179", Tm_179},
        {"Yb-179", Yb_179},
        {"Lu-179", Lu_179},
        {"Hf-179", Hf_179},
        {"Ta-179", Ta_179},
        {"W-179", W_179},
        {"Re-179", Re_179},
        {"Os-179", Os_179},
        {"Ir-179", Ir_179},
        {"Pt-179", Pt_179},
        {"Au-179", Au_179},
        {"Hg-179", Hg_179},
        {"Tl-179", Tl_179},
        {"Pb-179", Pb_179},
        {"Er-180", Er_180},
        {"Tm-180", Tm_180},
        {"Yb-180", Yb_180},
        {"Lu-180", Lu_180},
        {"Hf-180", Hf_180},
        {"Ta-180", Ta_180},
        {"W-180", W_180},
        {"Re-180", Re_180},
        {"Os-180", Os_180},
        {"Ir-180", Ir_180},
        {"Pt-180", Pt_180},
        {"Au-180", Au_180},
        {"Hg-180", Hg_180},
        {"Tl-180", Tl_180},
        {"Pb-180", Pb_180},
        {"Tm-181", Tm_181},
        {"Yb-181", Yb_181},
        {"Lu-181", Lu_181},
        {"Hf-181", Hf_181},
        {"Ta-181", Ta_181},
        {"W-181", W_181},
        {"Re-181", Re_181},
        {"Os-181", Os_181},
        {"Ir-181", Ir_181},
        {"Pt-181", Pt_181},
        {"Au-181", Au_181},
        {"Hg-181", Hg_181},
        {"Tl-181", Tl_181},
        {"Pb-181", Pb_181},
        {"Tm-182", Tm_182},
        {"Yb-182", Yb_182},
        {"Lu-182", Lu_182},
        {"Hf-182", Hf_182},
        {"Ta-182", Ta_182},
        {"W-182", W_182},
        {"Re-182", Re_182},
        {"Os-182", Os_182},
        {"Ir-182", Ir_182},
        {"Pt-182", Pt_182},
        {"Au-182", Au_182},
        {"Hg-182", Hg_182},
        {"Tl-182", Tl_182},
        {"Pb-182", Pb_182},
        {"Yb-183", Yb_183},
        {"Lu-183", Lu_183},
        {"Hf-183", Hf_183},
        {"Ta-183", Ta_183},
        {"W-183", W_183},
        {"Re-183", Re_183},
        {"Os-183", Os_183},
        {"Ir-183", Ir_183},
        {"Pt-183", Pt_183},
        {"Au-183", Au_183},
        {"Hg-183", Hg_183},
        {"Tl-183", Tl_183},
        {"Pb-183", Pb_183},
        {"Yb-184", Yb_184},
        {"Lu-184", Lu_184},
        {"Hf-184", Hf_184},
        {"Ta-184", Ta_184},
        {"W-184", W_184},
        {"Re-184", Re_184},
        {"Os-184", Os_184},
        {"Ir-184", Ir_184},
        {"Pt-184", Pt_184},
        {"Au-184", Au_184},
        {"Hg-184", Hg_184},
        {"Tl-184", Tl_184},
        {"Pb-184", Pb_184},
        {"Bi-184", Bi_184},
        {"Yb-185", Yb_185},
        {"Lu-185", Lu_185},
        {"Hf-185", Hf_185},
        {"Ta-185", Ta_185},
        {"W-185", W_185},
        {"Re-185", Re_185},
        {"Os-185", Os_185},
        {"Ir-185", Ir_185},
        {"Pt-185", Pt_185},
        {"Au-185", Au_185},
        {"Hg-185", Hg_185},
        {"Tl-185", Tl_185},
        {"Pb-185", Pb_185},
        {"Bi-185", Bi_185},
        {"Lu-186", Lu_186},
        {"Hf-186", Hf_186},
        {"Ta-186", Ta_186},
        {"W-186", W_186},
        {"Re-186", Re_186},
        {"Os-186", Os_186},
        {"Ir-186", Ir_186},
        {"Pt-186", Pt_186},
        {"Au-186", Au_186},
        {"Hg-186", Hg_186},
        {"Tl-186", Tl_186},
        {"Pb-186", Pb_186},
        {"Bi-186", Bi_186},
        {"Po-186", Po_186},
        {"Lu-187", Lu_187},
        {"Hf-187", Hf_187},
        {"Ta-187", Ta_187},
        {"W-187", W_187},
        {"Re-187", Re_187},
        {"Os-187", Os_187},
        {"Ir-187", Ir_187},
        {"Pt-187", Pt_187},
        {"Au-187", Au_187},
        {"Hg-187", Hg_187},
        {"Tl-187", Tl_187},
        {"Pb-187", Pb_187},
        {"Bi-187", Bi_187},
        {"Po-187", Po_187},
        {"Lu-188", Lu_188},
        {"Hf-188", Hf_188},
        {"Ta-188", Ta_188},
        {"W-188", W_188},
        {"Re-188", Re_188},
        {"Os-188", Os_188},
        {"Ir-188", Ir_188},
        {"Pt-188", Pt_188},
        {"Au-188", Au_188},
        {"Hg-188", Hg_188},
        {"Tl-188", Tl_188},
        {"Pb-188", Pb_188},
        {"Bi-188", Bi_188},
        {"Po-188", Po_188},
        {"Hf-189", Hf_189},
        {"Ta-189", Ta_189},
        {"W-189", W_189},
        {"Re-189", Re_189},
        {"Os-189", Os_189},
        {"Ir-189", Ir_189},
        {"Pt-189", Pt_189},
        {"Au-189", Au_189},
        {"Hg-189", Hg_189},
        {"Tl-189", Tl_189},
        {"Pb-189", Pb_189},
        {"Bi-189", Bi_189},
        {"Po-189", Po_189},
        {"Hf-190", Hf_190},
        {"Ta-190", Ta_190},
        {"W-190", W_190},
        {"Re-190", Re_190},
        {"Os-190", Os_190},
        {"Ir-190", Ir_190},
        {"Pt-190", Pt_190},
        {"Au-190", Au_190},
        {"Hg-190", Hg_190},
        {"Tl-190", Tl_190},
        {"Pb-190", Pb_190},
        {"Bi-190", Bi_190},
        {"Po-190", Po_190},
        {"Ta-191", Ta_191},
        {"W-191", W_191},
        {"Re-191", Re_191},
        {"Os-191", Os_191},
        {"Ir-191", Ir_191},
        {"Pt-191", Pt_191},
        {"Au-191", Au_191},
        {"Hg-191", Hg_191},
        {"Tl-191", Tl_191},
        {"Pb-191", Pb_191},
        {"Bi-191", Bi_191},
        {"Po-191", Po_191},
        {"At-191", At_191},
        {"Ta-192", Ta_192},
        {"W-192", W_192},
        {"Re-192", Re_192},
        {"Os-192", Os_192},
        {"Ir-192", Ir_192},
        {"Pt-192", Pt_192},
        {"Au-192", Au_192},
        {"Hg-192", Hg_192},
        {"Tl-192", Tl_192},
        {"Pb-192", Pb_192},
        {"Bi-192", Bi_192},
        {"Po-192", Po_192},
        {"At-192", At_192},
        {"Ta-193", Ta_193},
        {"W-193", W_193},
        {"Re-193", Re_193},
        {"Os-193", Os_193},
        {"Ir-193", Ir_193},
        {"Pt-193", Pt_193},
        {"Au-193", Au_193},
        {"Hg-193", Hg_193},
        {"Tl-193", Tl_193},
        {"Pb-193", Pb_193},
        {"Bi-193", Bi_193},
        {"Po-193", Po_193},
        {"At-193", At_193},
        {"Rn-193", Rn_193},
        {"Ta-194", Ta_194},
        {"W-194", W_194},
        {"Re-194", Re_194},
        {"Os-194", Os_194},
        {"Ir-194", Ir_194},
        {"Pt-194", Pt_194},
        {"Au-194", Au_194},
        {"Hg-194", Hg_194},
        {"Tl-194", Tl_194},
        {"Pb-194", Pb_194},
        {"Bi-194", Bi_194},
        {"Po-194", Po_194},
        {"At-194", At_194},
        {"Rn-194", Rn_194},
        {"W-195", W_195},
        {"Re-195", Re_195},
        {"Os-195", Os_195},
        {"Ir-195", Ir_195},
        {"Pt-195", Pt_195},
        {"Au-195", Au_195},
        {"Hg-195", Hg_195},
        {"Tl-195", Tl_195},
        {"Pb-195", Pb_195},
        {"Bi-195", Bi_195},
        {"Po-195", Po_195},
        {"At-195", At_195},
        {"Rn-195", Rn_195},
        {"W-196", W_196},
        {"Re-196", Re_196},
        {"Os-196", Os_196},
        {"Ir-196", Ir_196},
        {"Pt-196", Pt_196},
        {"Au-196", Au_196},
        {"Hg-196", Hg_196},
        {"Tl-196", Tl_196},
        {"Pb-196", Pb_196},
        {"Bi-196", Bi_196},
        {"Po-196", Po_196},
        {"At-196", At_196},
        {"Rn-196", Rn_196},
        {"W-197", W_197},
        {"Re-197", Re_197},
        {"Os-197", Os_197},
        {"Ir-197", Ir_197},
        {"Pt-197", Pt_197},
        {"Au-197", Au_197},
        {"Hg-197", Hg_197},
        {"Tl-197", Tl_197},
        {"Pb-197", Pb_197},
        {"Bi-197", Bi_197},
        {"Po-197", Po_197},
        {"At-197", At_197},
        {"Rn-197", Rn_197},
        {"Fr-197", Fr_197},
        {"Re-198", Re_198},
        {"Os-198", Os_198},
        {"Ir-198", Ir_198},
        {"Pt-198", Pt_198},
        {"Au-198", Au_198},
        {"Hg-198", Hg_198},
        {"Tl-198", Tl_198},
        {"Pb-198", Pb_198},
        {"Bi-198", Bi_198},
        {"Po-198", Po_198},
        {"At-198", At_198},
        {"Rn-198", Rn_198},
        {"Fr-198", Fr_198},
        {"Re-199", Re_199},
        {"Os-199", Os_199},
        {"Ir-199", Ir_199},
        {"Pt-199", Pt_199},
        {"Au-199", Au_199},
        {"Hg-199", Hg_199},
        {"Tl-199", Tl_199},
        {"Pb-199", Pb_199},
        {"Bi-199", Bi_199},
        {"Po-199", Po_199},
        {"At-199", At_199},
        {"Rn-199", Rn_199},
        {"Fr-199", Fr_199},
        {"Os-200", Os_200},
        {"Ir-200", Ir_200},
        {"Pt-200", Pt_200},
        {"Au-200", Au_200},
        {"Hg-200", Hg_200},
        {"Tl-200", Tl_200},
        {"Pb-200", Pb_200},
        {"Bi-200", Bi_200},
        {"Po-200", Po_200},
        {"At-200", At_200},
        {"Rn-200", Rn_200},
        {"Fr-200", Fr_200},
        {"Os-201", Os_201},
        {"Ir-201", Ir_201},
        {"Pt-201", Pt_201},
        {"Au-201", Au_201},
        {"Hg-201", Hg_201},
        {"Tl-201", Tl_201},
        {"Pb-201", Pb_201},
        {"Bi-201", Bi_201},
        {"Po-201", Po_201},
        {"At-201", At_201},
        {"Rn-201", Rn_201},
        {"Fr-201", Fr_201},
        {"Ra-201", Ra_201},
        {"Os-202", Os_202},
        {"Ir-202", Ir_202},
        {"Pt-202", Pt_202},
        {"Au-202", Au_202},
        {"Hg-202", Hg_202},
        {"Tl-202", Tl_202},
        {"Pb-202", Pb_202},
        {"Bi-202", Bi_202},
        {"Po-202", Po_202},
        {"At-202", At_202},
        {"Rn-202", Rn_202},
        {"Fr-202", Fr_202},
        {"Ra-202", Ra_202},
        {"Os-203", Os_203},
        {"Ir-203", Ir_203},
        {"Pt-203", Pt_203},
        {"Au-203", Au_203},
        {"Hg-203", Hg_203},
        {"Tl-203", Tl_203},
        {"Pb-203", Pb_203},
        {"Bi-203", Bi_203},
        {"Po-203", Po_203},
        {"At-203", At_203},
        {"Rn-203", Rn_203},
        {"Fr-203", Fr_203},
        {"Ra-203", Ra_203},
        {"Ir-204", Ir_204},
        {"Pt-204", Pt_204},
        {"Au-204", Au_204},
        {"Hg-204", Hg_204},
        {"Tl-204", Tl_204},
        {"Pb-204", Pb_204},
        {"Bi-204", Bi_204},
        {"Po-204", Po_204},
        {"At-204", At_204},
        {"Rn-204", Rn_204},
        {"Fr-204", Fr_204},
        {"Ra-204", Ra_204},
        {"Ir-205", Ir_205},
        {"Pt-205", Pt_205},
        {"Au-205", Au_205},
        {"Hg-205", Hg_205},
        {"Tl-205", Tl_205},
        {"Pb-205", Pb_205},
        {"Bi-205", Bi_205},
        {"Po-205", Po_205},
        {"At-205", At_205},
        {"Rn-205", Rn_205},
        {"Fr-205", Fr_205},
        {"Ra-205", Ra_205},
        {"Ac-205", Ac_205},
        {"Pt-206", Pt_206},
        {"Au-206", Au_206},
        {"Hg-206", Hg_206},
        {"Tl-206", Tl_206},
        {"Pb-206", Pb_206},
        {"Bi-206", Bi_206},
        {"Po-206", Po_206},
        {"At-206", At_206},
        {"Rn-206", Rn_206},
        {"Fr-206", Fr_206},
        {"Ra-206", Ra_206},
        {"Ac-206", Ac_206},
        {"Pt-207", Pt_207},
        {"Au-207", Au_207},
        {"Hg-207", Hg_207},
        {"Tl-207", Tl_207},
        {"Pb-207", Pb_207},
        {"Bi-207", Bi_207},
        {"Po-207", Po_207},
        {"At-207", At_207},
        {"Rn-207", Rn_207},
        {"Fr-207", Fr_207},
        {"Ra-207", Ra_207},
        {"Ac-207", Ac_207},
        {"Pt-208", Pt_208},
        {"Au-208", Au_208},
        {"Hg-208", Hg_208},
        {"Tl-208", Tl_208},
        {"Pb-208", Pb_208},
        {"Bi-208", Bi_208},
        {"Po-208", Po_208},
        {"At-208", At_208},
        {"Rn-208", Rn_208},
        {"Fr-208", Fr_208},
        {"Ra-208", Ra_208},
        {"Ac-208", Ac_208},
        {"Th-208", Th_208},
        {"Au-209", Au_209},
        {"Hg-209", Hg_209},
        {"Tl-209", Tl_209},
        {"Pb-209", Pb_209},
        {"Bi-209", Bi_209},
        {"Po-209", Po_209},
        {"At-209", At_209},
        {"Rn-209", Rn_209},
        {"Fr-209", Fr_209},
        {"Ra-209", Ra_209},
        {"Ac-209", Ac_209},
        {"Th-209", Th_209},
        {"Au-210", Au_210},
        {"Hg-210", Hg_210},
        {"Tl-210", Tl_210},
        {"Pb-210", Pb_210},
        {"Bi-210", Bi_210},
        {"Po-210", Po_210},
        {"At-210", At_210},
        {"Rn-210", Rn_210},
        {"Fr-210", Fr_210},
        {"Ra-210", Ra_210},
        {"Ac-210", Ac_210},
        {"Th-210", Th_210},
        {"Hg-211", Hg_211},
        {"Tl-211", Tl_211},
        {"Pb-211", Pb_211},
        {"Bi-211", Bi_211},
        {"Po-211", Po_211},
        {"At-211", At_211},
        {"Rn-211", Rn_211},
        {"Fr-211", Fr_211},
        {"Ra-211", Ra_211},
        {"Ac-211", Ac_211},
        {"Th-211", Th_211},
        {"Pa-211", Pa_211},
        {"Hg-212", Hg_212},
        {"Tl-212", Tl_212},
        {"Pb-212", Pb_212},
        {"Bi-212", Bi_212},
        {"Po-212", Po_212},
        {"At-212", At_212},
        {"Rn-212", Rn_212},
        {"Fr-212", Fr_212},
        {"Ra-212", Ra_212},
        {"Ac-212", Ac_212},
        {"Th-212", Th_212},
        {"Pa-212", Pa_212},
        {"Hg-213", Hg_213},
        {"Tl-213", Tl_213},
        {"Pb-213", Pb_213},
        {"Bi-213", Bi_213},
        {"Po-213", Po_213},
        {"At-213", At_213},
        {"Rn-213", Rn_213},
        {"Fr-213", Fr_213},
        {"Ra-213", Ra_213},
        {"Ac-213", Ac_213},
        {"Th-213", Th_213},
        {"Pa-213", Pa_213},
        {"Hg-214", Hg_214},
        {"Tl-214", Tl_214},
        {"Pb-214", Pb_214},
        {"Bi-214", Bi_214},
        {"Po-214", Po_214},
        {"At-214", At_214},
        {"Rn-214", Rn_214},
        {"Fr-214", Fr_214},
        {"Ra-214", Ra_214},
        {"Ac-214", Ac_214},
        {"Th-214", Th_214},
        {"Pa-214", Pa_214},
        {"Hg-215", Hg_215},
        {"Tl-215", Tl_215},
        {"Pb-215", Pb_215},
        {"Bi-215", Bi_215},
        {"Po-215", Po_215},
        {"At-215", At_215},
        {"Rn-215", Rn_215},
        {"Fr-215", Fr_215},
        {"Ra-215", Ra_215},
        {"Ac-215", Ac_215},
        {"Th-215", Th_215},
        {"Pa-215", Pa_215},
        {"U-215", U_215},
        {"Hg-216", Hg_216},
        {"Tl-216", Tl_216},
        {"Pb-216", Pb_216},
        {"Bi-216", Bi_216},
        {"Po-216", Po_216},
        {"At-216", At_216},
        {"Rn-216", Rn_216},
        {"Fr-216", Fr_216},
        {"Ra-216", Ra_216},
        {"Ac-216", Ac_216},
        {"Th-216", Th_216},
        {"Pa-216", Pa_216},
        {"U-216", U_216},
        {"Tl-217", Tl_217},
        {"Pb-217", Pb_217},
        {"Bi-217", Bi_217},
        {"Po-217", Po_217},
        {"At-217", At_217},
        {"Rn-217", Rn_217},
        {"Fr-217", Fr_217},
        {"Ra-217", Ra_217},
        {"Ac-217", Ac_217},
        {"Th-217", Th_217},
        {"Pa-217", Pa_217},
        {"U-217", U_217},
        {"Tl-218", Tl_218},
        {"Pb-218", Pb_218},
        {"Bi-218", Bi_218},
        {"Po-218", Po_218},
        {"At-218", At_218},
        {"Rn-218", Rn_218},
        {"Fr-218", Fr_218},
        {"Ra-218", Ra_218},
        {"Ac-218", Ac_218},
        {"Th-218", Th_218},
        {"Pa-218", Pa_218},
        {"U-218", U_218},
        {"Pb-219", Pb_219},
        {"Bi-219", Bi_219},
        {"Po-219", Po_219},
        {"At-219", At_219},
        {"Rn-219", Rn_219},
        {"Fr-219", Fr_219},
        {"Ra-219", Ra_219},
        {"Ac-219", Ac_219},
        {"Th-219", Th_219},
        {"Pa-219", Pa_219},
        {"U-219", U_219},
        {"Np-219", Np_219},
        {"Pb-220", Pb_220},
        {"Bi-220", Bi_220},
        {"Po-220", Po_220},
        {"At-220", At_220},
        {"Rn-220", Rn_220},
        {"Fr-220", Fr_220},
        {"Ra-220", Ra_220},
        {"Ac-220", Ac_220},
        {"Th-220", Th_220},
        {"Pa-220", Pa_220},
        {"U-220", U_220},
        {"Np-220", Np_220},
        {"Bi-221", Bi_221},
        {"Po-221", Po_221},
        {"At-221", At_221},
        {"Rn-221", Rn_221},
        {"Fr-221", Fr_221},
        {"Ra-221", Ra_221},
        {"Ac-221", Ac_221},
        {"Th-221", Th_221},
        {"Pa-221", Pa_221},
        {"U-221", U_221},
        {"Np-221", Np_221},
        {"Pu-221", Pu_221},
        {"Bi-222", Bi_222},
        {"Po-222", Po_222},
        {"At-222", At_222},
        {"Rn-222", Rn_222},
        {"Fr-222", Fr_222},
        {"Ra-222", Ra_222},
        {"Ac-222", Ac_222},
        {"Th-222", Th_222},
        {"Pa-222", Pa_222},
        {"U-222", U_222},
        {"Np-222", Np_222},
        {"Pu-222", Pu_222},
        {"Bi-223", Bi_223},
        {"Po-223", Po_223},
        {"At-223", At_223},
        {"Rn-223", Rn_223},
        {"Fr-223", Fr_223},
        {"Ra-223", Ra_223},
        {"Ac-223", Ac_223},
        {"Th-223", Th_223},
        {"Pa-223", Pa_223},
        {"U-223", U_223},
        {"Np-223", Np_223},
        {"Pu-223", Pu_223},
        {"Am-223", Am_223},
        {"Bi-224", Bi_224},
        {"Po-224", Po_224},
        {"At-224", At_224},
        {"Rn-224", Rn_224},
        {"Fr-224", Fr_224},
        {"Ra-224", Ra_224},
        {"Ac-224", Ac_224},
        {"Th-224", Th_224},
        {"Pa-224", Pa_224},
        {"U-224", U_224},
        {"Np-224", Np_224},
        {"Pu-224", Pu_224},
        {"Am-224", Am_224},
        {"Po-225", Po_225},
        {"At-225", At_225},
        {"Rn-225", Rn_225},
        {"Fr-225", Fr_225},
        {"Ra-225", Ra_225},
        {"Ac-225", Ac_225},
        {"Th-225", Th_225},
        {"Pa-225", Pa_225},
        {"U-225", U_225},
        {"Np-225", Np_225},
        {"Pu-225", Pu_225},
        {"Am-225", Am_225},
        {"Po-226", Po_226},
        {"At-226", At_226},
        {"Rn-226", Rn_226},
        {"Fr-226", Fr_226},
        {"Ra-226", Ra_226},
        {"Ac-226", Ac_226},
        {"Th-226", Th_226},
        {"Pa-226", Pa_226},
        {"U-226", U_226},
        {"Np-226", Np_226},
        {"Pu-226", Pu_226},
        {"Am-226", Am_226},
        {"Po-227", Po_227},
        {"At-227", At_227},
        {"Rn-227", Rn_227},
        {"Fr-227", Fr_227},
        {"Ra-227", Ra_227},
        {"Ac-227", Ac_227},
        {"Th-227", Th_227},
        {"Pa-227", Pa_227},
        {"U-227", U_227},
        {"Np-227", Np_227},
        {"Pu-227", Pu_227},
        {"Am-227", Am_227},
        {"At-228", At_228},
        {"Rn-228", Rn_228},
        {"Fr-228", Fr_228},
        {"Ra-228", Ra_228},
        {"Ac-228", Ac_228},
        {"Th-228", Th_228},
        {"Pa-228", Pa_228},
        {"U-228", U_228},
        {"Np-228", Np_228},
        {"Pu-228", Pu_228},
        {"Am-228", Am_228},
        {"At-229", At_229},
        {"Rn-229", Rn_229},
        {"Fr-229", Fr_229},
        {"Ra-229", Ra_229},
        {"Ac-229", Ac_229},
        {"Th-229", Th_229},
        {"Pa-229", Pa_229},
        {"U-229", U_229},
        {"Np-229", Np_229},
        {"Pu-229", Pu_229},
        {"Am-229", Am_229},
        {"Rn-230", Rn_230},
        {"Fr-230", Fr_230},
        {"Ra-230", Ra_230},
        {"Ac-230", Ac_230},
        {"Th-230", Th_230},
        {"Pa-230", Pa_230},
        {"U-230", U_230},
        {"Np-230", Np_230},
        {"Pu-230", Pu_230},
        {"Am-230", Am_230},
        {"Rn-231", Rn_231},
        {"Fr-231", Fr_231},
        {"Ra-231", Ra_231},
        {"Ac-231", Ac_231},
        {"Th-231", Th_231},
        {"Pa-231", Pa_231},
        {"U-231", U_231},
        {"Np-231", Np_231},
        {"Pu-231", Pu_231},
        {"Am-231", Am_231},
        {"Cm-231", Cm_231},
        {"Fr-232", Fr_232},
        {"Ra-232", Ra_232},
        {"Ac-232", Ac_232},
        {"Th-232", Th_232},
        {"Pa-232", Pa_232},
        {"U-232", U_232},
        {"Np-232", Np_232},
        {"Pu-232", Pu_232},
        {"Am-232", Am_232},
        {"Cm-232", Cm_232},
        {"Fr-233", Fr_233},
        {"Ra-233", Ra_233},
        {"Ac-233", Ac_233},
        {"Th-233", Th_233},
        {"Pa-233", Pa_233},
        {"U-233", U_233},
        {"Np-233", Np_233},
        {"Pu-233", Pu_233},
        {"Am-233", Am_233},
        {"Cm-233", Cm_233},
        {"Bk-233", Bk_233},
        {"Ra-234", Ra_234},
        {"Ac-234", Ac_234},
        {"Th-234", Th_234},
        {"Pa-234", Pa_234},
        {"U-234", U_234},
        {"Np-234", Np_234},
        {"Pu-234", Pu_234},
        {"Am-234", Am_234},
        {"Cm-234", Cm_234},
        {"Bk-234", Bk_234},
        {"Ra-235", Ra_235},
        {"Ac-235", Ac_235},
        {"Th-235", Th_235},
        {"Pa-235", Pa_235},
        {"U-235", U_235},
        {"Np-235", Np_235},
        {"Pu-235", Pu_235},
        {"Am-235", Am_235},
        {"Cm-235", Cm_235},
        {"Bk-235", Bk_235},
        {"Ac-236", Ac_236},
        {"Th-236", Th_236},
        {"Pa-236", Pa_236},
        {"U-236", U_236},
        {"Np-236", Np_236},
        {"Pu-236", Pu_236},
        {"Am-236", Am_236},
        {"Cm-236", Cm_236},
        {"Bk-236", Bk_236},
        {"Ac-237", Ac_237},
        {"Th-237", Th_237},
        {"Pa-237", Pa_237},
        {"U-237", U_237},
        {"Np-237", Np_237},
        {"Pu-237", Pu_237},
        {"Am-237", Am_237},
        {"Cm-237", Cm_237},
        {"Bk-237", Bk_237},
        {"Cf-237", Cf_237},
        {"Th-238", Th_238},
        {"Pa-238", Pa_238},
        {"U-238", U_238},
        {"Np-238", Np_238},
        {"Pu-238", Pu_238},
        {"Am-238", Am_238},
        {"Cm-238", Cm_238},
        {"Bk-238", Bk_238},
        {"Cf-238", Cf_238},
        {"Th-239", Th_239},
        {"Pa-239", Pa_239},
        {"U-239", U_239},
        {"Np-239", Np_239},
        {"Pu-239", Pu_239},
        {"Am-239", Am_239},
        {"Cm-239", Cm_239},
        {"Bk-239", Bk_239},
        {"Cf-239", Cf_239},
        {"Es-239", Es_239},
        {"Pa-240", Pa_240},
        {"U-240", U_240},
        {"Np-240", Np_240},
        {"Pu-240", Pu_240},
        {"Am-240", Am_240},
        {"Cm-240", Cm_240},
        {"Bk-240", Bk_240},
        {"Cf-240", Cf_240},
        {"Es-240", Es_240},
        {"Pa-241", Pa_241},
        {"U-241", U_241},
        {"Np-241", Np_241},
        {"Pu-241", Pu_241},
        {"Am-241", Am_241},
        {"Cm-241", Cm_241},
        {"Bk-241", Bk_241},
        {"Cf-241", Cf_241},
        {"Es-241", Es_241},
        {"Fm-241", Fm_241},
        {"U-242", U_242},
        {"Np-242", Np_242},
        {"Pu-242", Pu_242},
        {"Am-242", Am_242},
        {"Cm-242", Cm_242},
        {"Bk-242", Bk_242},
        {"Cf-242", Cf_242},
        {"Es-242", Es_242},
        {"Fm-242", Fm_242},
        {"U-243", U_243},
        {"Np-243", Np_243},
        {"Pu-243", Pu_243},
        {"Am-243", Am_243},
        {"Cm-243", Cm_243},
        {"Bk-243", Bk_243},
        {"Cf-243", Cf_243},
        {"Es-243", Es_243},
        {"Fm-243", Fm_243},
        {"Np-244", Np_244},
        {"Pu-244", Pu_244},
        {"Am-244", Am_244},
        {"Cm-244", Cm_244},
        {"Bk-244", Bk_244},
        {"Cf-244", Cf_244},
        {"Es-244", Es_244},
        {"Fm-244", Fm_244},
        {"Md-244", Md_244},
        {"Np-245", Np_245},
        {"Pu-245", Pu_245},
        {"Am-245", Am_245},
        {"Cm-245", Cm_245},
        {"Bk-245", Bk_245},
        {"Cf-245", Cf_245},
        {"Es-245", Es_245},
        {"Fm-245", Fm_245},
        {"Md-245", Md_245},
        {"Pu-246", Pu_246},
        {"Am-246", Am_246},
        {"Cm-246", Cm_246},
        {"Bk-246", Bk_246},
        {"Cf-246", Cf_246},
        {"Es-246", Es_246},
        {"Fm-246", Fm_246},
        {"Md-246", Md_246},
        {"Pu-247", Pu_247},
        {"Am-247", Am_247},
        {"Cm-247", Cm_247},
        {"Bk-247", Bk_247},
        {"Cf-247", Cf_247},
        {"Es-247", Es_247},
        {"Fm-247", Fm_247},
        {"Md-247", Md_247},
        {"Am-248", Am_248},
        {"Cm-248", Cm_248},
        {"Bk-248", Bk_248},
        {"Cf-248", Cf_248},
        {"Es-248", Es_248},
        {"Fm-248", Fm_248},
        {"Md-248", Md_248},
        {"No-248", No_248},
        {"Am-249", Am_249},
        {"Cm-249", Cm_249},
        {"Bk-249", Bk_249},
        {"Cf-249", Cf_249},
        {"Es-249", Es_249},
        {"Fm-249", Fm_249},
        {"Md-249", Md_249},
        {"No-249", No_249},
        {"Cm-250", Cm_250},
        {"Bk-250", Bk_250},
        {"Cf-250", Cf_250},
        {"Es-250", Es_250},
        {"Fm-250", Fm_250},
        {"Md-250", Md_250},
        {"No-250", No_250},
        {"Cm-251", Cm_251},
        {"Bk-251", Bk_251},
        {"Cf-251", Cf_251},
        {"Es-251", Es_251},
        {"Fm-251", Fm_251},
        {"Md-251", Md_251},
        {"No-251", No_251},
        {"Lr-251", Lr_251},
        {"Cm-252", Cm_252},
        {"Bk-252", Bk_252},
        {"Cf-252", Cf_252},
        {"Es-252", Es_252},
        {"Fm-252", Fm_252},
        {"Md-252", Md_252},
        {"No-252", No_252},
        {"Lr-252", Lr_252},
        {"Bk-253", Bk_253},
        {"Cf-253", Cf_253},
        {"Es-253", Es_253},
        {"Fm-253", Fm_253},
        {"Md-253", Md_253},
        {"No-253", No_253},
        {"Lr-253", Lr_253},
        {"Rf-253", Rf_253},
        {"Bk-254", Bk_254},
        {"Cf-254", Cf_254},
        {"Es-254", Es_254},
        {"Fm-254", Fm_254},
        {"Md-254", Md_254},
        {"No-254", No_254},
        {"Lr-254", Lr_254},
        {"Rf-254", Rf_254},
        {"Cf-255", Cf_255},
        {"Es-255", Es_255},
        {"Fm-255", Fm_255},
        {"Md-255", Md_255},
        {"No-255", No_255},
        {"Lr-255", Lr_255},
        {"Rf-255", Rf_255},
        {"Db-255", Db_255},
        {"Cf-256", Cf_256},
        {"Es-256", Es_256},
        {"Fm-256", Fm_256},
        {"Md-256", Md_256},
        {"No-256", No_256},
        {"Lr-256", Lr_256},
        {"Rf-256", Rf_256},
        {"Db-256", Db_256},
        {"Es-257", Es_257},
        {"Fm-257", Fm_257},
        {"Md-257", Md_257},
        {"No-257", No_257},
        {"Lr-257", Lr_257},
        {"Rf-257", Rf_257},
        {"Db-257", Db_257},
        {"Es-258", Es_258},
        {"Fm-258", Fm_258},
        {"Md-258", Md_258},
        {"No-258", No_258},
        {"Lr-258", Lr_258},
        {"Rf-258", Rf_258},
        {"Db-258", Db_258},
        {"Sg-258", Sg_258},
        {"Fm-259", Fm_259},
        {"Md-259", Md_259},
        {"No-259", No_259},
        {"Lr-259", Lr_259},
        {"Rf-259", Rf_259},
        {"Db-259", Db_259},
        {"Sg-259", Sg_259},
        {"Fm-260", Fm_260},
        {"Md-260", Md_260},
        {"No-260", No_260},
        {"Lr-260", Lr_260},
        {"Rf-260", Rf_260},
        {"Db-260", Db_260},
        {"Sg-260", Sg_260},
        {"Bh-260", Bh_260},
        {"Md-261", Md_261},
        {"No-261", No_261},
        {"Lr-261", Lr_261},
        {"Rf-261", Rf_261},
        {"Db-261", Db_261},
        {"Sg-261", Sg_261},
        {"Bh-261", Bh_261},
        {"Md-262", Md_262},
        {"No-262", No_262},
        {"Lr-262", Lr_262},
        {"Rf-262", Rf_262},
        {"Db-262", Db_262},
        {"Sg-262", Sg_262},
        {"Bh-262", Bh_262},
        {"No-263", No_263},
        {"Lr-263", Lr_263},
        {"Rf-263", Rf_263},
        {"Db-263", Db_263},
        {"Sg-263", Sg_263},
        {"Bh-263", Bh_263},
        {"Hs-263", Hs_263},
        {"No-264", No_264},
        {"Lr-264", Lr_264},
        {"Rf-264", Rf_264},
        {"Db-264", Db_264},
        {"Sg-264", Sg_264},
        {"Bh-264", Bh_264},
        {"Hs-264", Hs_264},
        {"Lr-265", Lr_265},
        {"Rf-265", Rf_265},
        {"Db-265", Db_265},
        {"Sg-265", Sg_265},
        {"Bh-265", Bh_265},
        {"Hs-265", Hs_265},
        {"Mt-265", Mt_265},
        {"Lr-266", Lr_266},
        {"Rf-266", Rf_266},
        {"Db-266", Db_266},
        {"Sg-266", Sg_266},
        {"Bh-266", Bh_266},
        {"Hs-266", Hs_266},
        {"Mt-266", Mt_266},
        {"Rf-267", Rf_267},
        {"Db-267", Db_267},
        {"Sg-267", Sg_267},
        {"Bh-267", Bh_267},
        {"Hs-267", Hs_267},
        {"Mt-267", Mt_267},
        {"Ds-267", Ds_267},
        {"Rf-268", Rf_268},
        {"Db-268", Db_268},
        {"Sg-268", Sg_268},
        {"Bh-268", Bh_268},
        {"Hs-268", Hs_268},
        {"Mt-268", Mt_268},
        {"Ds-268", Ds_268},
        {"Db-269", Db_269},
        {"Sg-269", Sg_269},
        {"Bh-269", Bh_269},
        {"Hs-269", Hs_269},
        {"Mt-269", Mt_269},
        {"Ds-269", Ds_269},
        {"Db-270", Db_270},
        {"Sg-270", Sg_270},
        {"Bh-270", Bh_270},
        {"Hs-270", Hs_270},
        {"Mt-270", Mt_270},
        {"Ds-270", Ds_270},
        {"Sg-271", Sg_271},
        {"Bh-271", Bh_271},
        {"Hs-271", Hs_271},
        {"Mt-271", Mt_271},
        {"Ds-271", Ds_271},
        {"Sg-272", Sg_272},
        {"Bh-272", Bh_272},
        {"Hs-272", Hs_272},
        {"Mt-272", Mt_272},
        {"Ds-272", Ds_272},
        {"Rg-272", Rg_272},
        {"Sg-273", Sg_273},
        {"Bh-273", Bh_273},
        {"Hs-273", Hs_273},
        {"Mt-273", Mt_273},
        {"Ds-273", Ds_273},
        {"Rg-273", Rg_273},
        {"Bh-274", Bh_274},
        {"Hs-274", Hs_274},
        {"Mt-274", Mt_274},
        {"Ds-274", Ds_274},
        {"Rg-274", Rg_274},
        {"Bh-275", Bh_275},
        {"Hs-275", Hs_275},
        {"Mt-275", Mt_275},
        {"Ds-275", Ds_275},
        {"Rg-275", Rg_275},
        {"Bh-276", Bh_276},
        {"Hs-276", Hs_276},
        {"Mt-276", Mt_276},
        {"Ds-276", Ds_276},
        {"Rg-276", Rg_276},
        {"Cn-276", Cn_276},
        {"Bh-277", Bh_277},
        {"Hs-277", Hs_277},
        {"Mt-277", Mt_277},
        {"Ds-277", Ds_277},
        {"Rg-277", Rg_277},
        {"Cn-277", Cn_277},
        {"Bh-278", Bh_278},
        {"Hs-278", Hs_278},
        {"Mt-278", Mt_278},
        {"Ds-278", Ds_278},
        {"Rg-278", Rg_278},
        {"Cn-278", Cn_278},
        {"Nh-278", Nh_278},
        {"Hs-279", Hs_279},
        {"Mt-279", Mt_279},
        {"Ds-279", Ds_279},
        {"Rg-279", Rg_279},
        {"Cn-279", Cn_279},
        {"Nh-279", Nh_279},
        {"Hs-280", Hs_280},
        {"Mt-280", Mt_280},
        {"Ds-280", Ds_280},
        {"Rg-280", Rg_280},
        {"Cn-280", Cn_280},
        {"Nh-280", Nh_280},
        {"Mt-281", Mt_281},
        {"Ds-281", Ds_281},
        {"Rg-281", Rg_281},
        {"Cn-281", Cn_281},
        {"Nh-281", Nh_281},
        {"Mt-282", Mt_282},
        {"Ds-282", Ds_282},
        {"Rg-282", Rg_282},
        {"Cn-282", Cn_282},
        {"Nh-282", Nh_282},
        {"Ds-283", Ds_283},
        {"Rg-283", Rg_283},
        {"Cn-283", Cn_283},
        {"Nh-283", Nh_283},
        {"Ds-284", Ds_284},
        {"Rg-284", Rg_284},
        {"Cn-284", Cn_284},
        {"Nh-284", Nh_284},
        {"Fl-284", Fl_284},
        {"Rg-285", Rg_285},
        {"Cn-285", Cn_285},
        {"Nh-285", Nh_285},
        {"Fl-285", Fl_285},
        {"Rg-286", Rg_286},
        {"Cn-286", Cn_286},
        {"Nh-286", Nh_286},
        {"Fl-286", Fl_286},
        {"Cn-287", Cn_287},
        {"Nh-287", Nh_287},
        {"Fl-287", Fl_287},
        {"Mc-287", Mc_287},
        {"Cn-288", Cn_288},
        {"Nh-288", Nh_288},
        {"Fl-288", Fl_288},
        {"Mc-288", Mc_288},
        {"Nh-289", Nh_289},
        {"Fl-289", Fl_289},
        {"Mc-289", Mc_289},
        {"Lv-289", Lv_289},
        {"Nh-290", Nh_290},
        {"Fl-290", Fl_290},
        {"Mc-290", Mc_290},
        {"Lv-290", Lv_290},
        {"Fl-291", Fl_291},
        {"Mc-291", Mc_291},
        {"Lv-291", Lv_291},
        {"Ts-291", Ts_291},
        {"Mc-292", Mc_292},
        {"Lv-292", Lv_292},
        {"Ts-292", Ts_292},
        {"Lv-293", Lv_293},
        {"Ts-293", Ts_293},
        {"Og-293", Og_293},
        {"Ts-294", Ts_294},
        {"Og-294", Og_294},
        {"Og-295", Og_295},
    };
}; // namespace chemSpecies

#ifndef N_1
    #define N_1
#endif // N_1
#ifndef H_1
    #define H_1
#endif // H_1
#ifndef H_2
    #define H_2
#endif // H_2
#ifndef H_3
    #define H_3
#endif // H_3
#ifndef HE_3
    #define HE_3
#endif // HE_3
#ifndef LI_3
    #define LI_3
#endif // LI_3
#ifndef H_4
    #define H_4
#endif // H_4
#ifndef HE_4
    #define HE_4
#endif // HE_4
#ifndef LI_4
    #define LI_4
#endif // LI_4
#ifndef H_5
    #define H_5
#endif // H_5
#ifndef HE_5
    #define HE_5
#endif // HE_5
#ifndef LI_5
    #define LI_5
#endif // LI_5
#ifndef BE_5
    #define BE_5
#endif // BE_5
#ifndef H_6
    #define H_6
#endif // H_6
#ifndef HE_6
    #define HE_6
#endif // HE_6
#ifndef LI_6
    #define LI_6
#endif // LI_6
#ifndef BE_6
    #define BE_6
#endif // BE_6
#ifndef B_6
    #define B_6
#endif // B_6
#ifndef H_7
    #define H_7
#endif // H_7
#ifndef HE_7
    #define HE_7
#endif // HE_7
#ifndef LI_7
    #define LI_7
#endif // LI_7
#ifndef BE_7
    #define BE_7
#endif // BE_7
#ifndef B_7
    #define B_7
#endif // B_7
#ifndef HE_8
    #define HE_8
#endif // HE_8
#ifndef LI_8
    #define LI_8
#endif // LI_8
#ifndef BE_8
    #define BE_8
#endif // BE_8
#ifndef B_8
    #define B_8
#endif // B_8
#ifndef C_8
    #define C_8
#endif // C_8
#ifndef HE_9
    #define HE_9
#endif // HE_9
#ifndef LI_9
    #define LI_9
#endif // LI_9
#ifndef BE_9
    #define BE_9
#endif // BE_9
#ifndef B_9
    #define B_9
#endif // B_9
#ifndef C_9
    #define C_9
#endif // C_9
#ifndef HE_10
    #define HE_10
#endif // HE_10
#ifndef LI_10
    #define LI_10
#endif // LI_10
#ifndef BE_10
    #define BE_10
#endif // BE_10
#ifndef B_10
    #define B_10
#endif // B_10
#ifndef C_10
    #define C_10
#endif // C_10
#ifndef N_10
    #define N_10
#endif // N_10
#ifndef LI_11
    #define LI_11
#endif // LI_11
#ifndef BE_11
    #define BE_11
#endif // BE_11
#ifndef B_11
    #define B_11
#endif // B_11
#ifndef C_11
    #define C_11
#endif // C_11
#ifndef N_11
    #define N_11
#endif // N_11
#ifndef O_11
    #define O_11
#endif // O_11
#ifndef LI_12
    #define LI_12
#endif // LI_12
#ifndef BE_12
    #define BE_12
#endif // BE_12
#ifndef B_12
    #define B_12
#endif // B_12
#ifndef C_12
    #define C_12
#endif // C_12
#ifndef N_12
    #define N_12
#endif // N_12
#ifndef O_12
    #define O_12
#endif // O_12
#ifndef LI_13
    #define LI_13
#endif // LI_13
#ifndef BE_13
    #define BE_13
#endif // BE_13
#ifndef B_13
    #define B_13
#endif // B_13
#ifndef C_13
    #define C_13
#endif // C_13
#ifndef N_13
    #define N_13
#endif // N_13
#ifndef O_13
    #define O_13
#endif // O_13
#ifndef F_13
    #define F_13
#endif // F_13
#ifndef BE_14
    #define BE_14
#endif // BE_14
#ifndef B_14
    #define B_14
#endif // B_14
#ifndef C_14
    #define C_14
#endif // C_14
#ifndef N_14
    #define N_14
#endif // N_14
#ifndef O_14
    #define O_14
#endif // O_14
#ifndef F_14
    #define F_14
#endif // F_14
#ifndef BE_15
    #define BE_15
#endif // BE_15
#ifndef B_15
    #define B_15
#endif // B_15
#ifndef C_15
    #define C_15
#endif // C_15
#ifndef N_15
    #define N_15
#endif // N_15
#ifndef O_15
    #define O_15
#endif // O_15
#ifndef F_15
    #define F_15
#endif // F_15
#ifndef NE_15
    #define NE_15
#endif // NE_15
#ifndef BE_16
    #define BE_16
#endif // BE_16
#ifndef B_16
    #define B_16
#endif // B_16
#ifndef C_16
    #define C_16
#endif // C_16
#ifndef N_16
    #define N_16
#endif // N_16
#ifndef O_16
    #define O_16
#endif // O_16
#ifndef F_16
    #define F_16
#endif // F_16
#ifndef NE_16
    #define NE_16
#endif // NE_16
#ifndef B_17
    #define B_17
#endif // B_17
#ifndef C_17
    #define C_17
#endif // C_17
#ifndef N_17
    #define N_17
#endif // N_17
#ifndef O_17
    #define O_17
#endif // O_17
#ifndef F_17
    #define F_17
#endif // F_17
#ifndef NE_17
    #define NE_17
#endif // NE_17
#ifndef NA_17
    #define NA_17
#endif // NA_17
#ifndef B_18
    #define B_18
#endif // B_18
#ifndef C_18
    #define C_18
#endif // C_18
#ifndef N_18
    #define N_18
#endif // N_18
#ifndef O_18
    #define O_18
#endif // O_18
#ifndef F_18
    #define F_18
#endif // F_18
#ifndef NE_18
    #define NE_18
#endif // NE_18
#ifndef NA_18
    #define NA_18
#endif // NA_18
#ifndef B_19
    #define B_19
#endif // B_19
#ifndef C_19
    #define C_19
#endif // C_19
#ifndef N_19
    #define N_19
#endif // N_19
#ifndef O_19
    #define O_19
#endif // O_19
#ifndef F_19
    #define F_19
#endif // F_19
#ifndef NE_19
    #define NE_19
#endif // NE_19
#ifndef NA_19
    #define NA_19
#endif // NA_19
#ifndef MG_19
    #define MG_19
#endif // MG_19
#ifndef B_20
    #define B_20
#endif // B_20
#ifndef C_20
    #define C_20
#endif // C_20
#ifndef N_20
    #define N_20
#endif // N_20
#ifndef O_20
    #define O_20
#endif // O_20
#ifndef F_20
    #define F_20
#endif // F_20
#ifndef NE_20
    #define NE_20
#endif // NE_20
#ifndef NA_20
    #define NA_20
#endif // NA_20
#ifndef MG_20
    #define MG_20
#endif // MG_20
#ifndef B_21
    #define B_21
#endif // B_21
#ifndef C_21
    #define C_21
#endif // C_21
#ifndef N_21
    #define N_21
#endif // N_21
#ifndef O_21
    #define O_21
#endif // O_21
#ifndef F_21
    #define F_21
#endif // F_21
#ifndef NE_21
    #define NE_21
#endif // NE_21
#ifndef NA_21
    #define NA_21
#endif // NA_21
#ifndef MG_21
    #define MG_21
#endif // MG_21
#ifndef AL_21
    #define AL_21
#endif // AL_21
#ifndef C_22
    #define C_22
#endif // C_22
#ifndef N_22
    #define N_22
#endif // N_22
#ifndef O_22
    #define O_22
#endif // O_22
#ifndef F_22
    #define F_22
#endif // F_22
#ifndef NE_22
    #define NE_22
#endif // NE_22
#ifndef NA_22
    #define NA_22
#endif // NA_22
#ifndef MG_22
    #define MG_22
#endif // MG_22
#ifndef AL_22
    #define AL_22
#endif // AL_22
#ifndef SI_22
    #define SI_22
#endif // SI_22
#ifndef C_23
    #define C_23
#endif // C_23
#ifndef N_23
    #define N_23
#endif // N_23
#ifndef O_23
    #define O_23
#endif // O_23
#ifndef F_23
    #define F_23
#endif // F_23
#ifndef NE_23
    #define NE_23
#endif // NE_23
#ifndef NA_23
    #define NA_23
#endif // NA_23
#ifndef MG_23
    #define MG_23
#endif // MG_23
#ifndef AL_23
    #define AL_23
#endif // AL_23
#ifndef SI_23
    #define SI_23
#endif // SI_23
#ifndef N_24
    #define N_24
#endif // N_24
#ifndef O_24
    #define O_24
#endif // O_24
#ifndef F_24
    #define F_24
#endif // F_24
#ifndef NE_24
    #define NE_24
#endif // NE_24
#ifndef NA_24
    #define NA_24
#endif // NA_24
#ifndef MG_24
    #define MG_24
#endif // MG_24
#ifndef AL_24
    #define AL_24
#endif // AL_24
#ifndef SI_24
    #define SI_24
#endif // SI_24
#ifndef P_24
    #define P_24
#endif // P_24
#ifndef N_25
    #define N_25
#endif // N_25
#ifndef O_25
    #define O_25
#endif // O_25
#ifndef F_25
    #define F_25
#endif // F_25
#ifndef NE_25
    #define NE_25
#endif // NE_25
#ifndef NA_25
    #define NA_25
#endif // NA_25
#ifndef MG_25
    #define MG_25
#endif // MG_25
#ifndef AL_25
    #define AL_25
#endif // AL_25
#ifndef SI_25
    #define SI_25
#endif // SI_25
#ifndef P_25
    #define P_25
#endif // P_25
#ifndef O_26
    #define O_26
#endif // O_26
#ifndef F_26
    #define F_26
#endif // F_26
#ifndef NE_26
    #define NE_26
#endif // NE_26
#ifndef NA_26
    #define NA_26
#endif // NA_26
#ifndef MG_26
    #define MG_26
#endif // MG_26
#ifndef AL_26
    #define AL_26
#endif // AL_26
#ifndef SI_26
    #define SI_26
#endif // SI_26
#ifndef P_26
    #define P_26
#endif // P_26
#ifndef S_26
    #define S_26
#endif // S_26
#ifndef O_27
    #define O_27
#endif // O_27
#ifndef F_27
    #define F_27
#endif // F_27
#ifndef NE_27
    #define NE_27
#endif // NE_27
#ifndef NA_27
    #define NA_27
#endif // NA_27
#ifndef MG_27
    #define MG_27
#endif // MG_27
#ifndef AL_27
    #define AL_27
#endif // AL_27
#ifndef SI_27
    #define SI_27
#endif // SI_27
#ifndef P_27
    #define P_27
#endif // P_27
#ifndef S_27
    #define S_27
#endif // S_27
#ifndef O_28
    #define O_28
#endif // O_28
#ifndef F_28
    #define F_28
#endif // F_28
#ifndef NE_28
    #define NE_28
#endif // NE_28
#ifndef NA_28
    #define NA_28
#endif // NA_28
#ifndef MG_28
    #define MG_28
#endif // MG_28
#ifndef AL_28
    #define AL_28
#endif // AL_28
#ifndef SI_28
    #define SI_28
#endif // SI_28
#ifndef P_28
    #define P_28
#endif // P_28
#ifndef S_28
    #define S_28
#endif // S_28
#ifndef CL_28
    #define CL_28
#endif // CL_28
#ifndef F_29
    #define F_29
#endif // F_29
#ifndef NE_29
    #define NE_29
#endif // NE_29
#ifndef NA_29
    #define NA_29
#endif // NA_29
#ifndef MG_29
    #define MG_29
#endif // MG_29
#ifndef AL_29
    #define AL_29
#endif // AL_29
#ifndef SI_29
    #define SI_29
#endif // SI_29
#ifndef P_29
    #define P_29
#endif // P_29
#ifndef S_29
    #define S_29
#endif // S_29
#ifndef CL_29
    #define CL_29
#endif // CL_29
#ifndef AR_29
    #define AR_29
#endif // AR_29
#ifndef F_30
    #define F_30
#endif // F_30
#ifndef NE_30
    #define NE_30
#endif // NE_30
#ifndef NA_30
    #define NA_30
#endif // NA_30
#ifndef MG_30
    #define MG_30
#endif // MG_30
#ifndef AL_30
    #define AL_30
#endif // AL_30
#ifndef SI_30
    #define SI_30
#endif // SI_30
#ifndef P_30
    #define P_30
#endif // P_30
#ifndef S_30
    #define S_30
#endif // S_30
#ifndef CL_30
    #define CL_30
#endif // CL_30
#ifndef AR_30
    #define AR_30
#endif // AR_30
#ifndef F_31
    #define F_31
#endif // F_31
#ifndef NE_31
    #define NE_31
#endif // NE_31
#ifndef NA_31
    #define NA_31
#endif // NA_31
#ifndef MG_31
    #define MG_31
#endif // MG_31
#ifndef AL_31
    #define AL_31
#endif // AL_31
#ifndef SI_31
    #define SI_31
#endif // SI_31
#ifndef P_31
    #define P_31
#endif // P_31
#ifndef S_31
    #define S_31
#endif // S_31
#ifndef CL_31
    #define CL_31
#endif // CL_31
#ifndef AR_31
    #define AR_31
#endif // AR_31
#ifndef K_31
    #define K_31
#endif // K_31
#ifndef NE_32
    #define NE_32
#endif // NE_32
#ifndef NA_32
    #define NA_32
#endif // NA_32
#ifndef MG_32
    #define MG_32
#endif // MG_32
#ifndef AL_32
    #define AL_32
#endif // AL_32
#ifndef SI_32
    #define SI_32
#endif // SI_32
#ifndef P_32
    #define P_32
#endif // P_32
#ifndef S_32
    #define S_32
#endif // S_32
#ifndef CL_32
    #define CL_32
#endif // CL_32
#ifndef AR_32
    #define AR_32
#endif // AR_32
#ifndef K_32
    #define K_32
#endif // K_32
#ifndef NE_33
    #define NE_33
#endif // NE_33
#ifndef NA_33
    #define NA_33
#endif // NA_33
#ifndef MG_33
    #define MG_33
#endif // MG_33
#ifndef AL_33
    #define AL_33
#endif // AL_33
#ifndef SI_33
    #define SI_33
#endif // SI_33
#ifndef P_33
    #define P_33
#endif // P_33
#ifndef S_33
    #define S_33
#endif // S_33
#ifndef CL_33
    #define CL_33
#endif // CL_33
#ifndef AR_33
    #define AR_33
#endif // AR_33
#ifndef K_33
    #define K_33
#endif // K_33
#ifndef CA_33
    #define CA_33
#endif // CA_33
#ifndef NE_34
    #define NE_34
#endif // NE_34
#ifndef NA_34
    #define NA_34
#endif // NA_34
#ifndef MG_34
    #define MG_34
#endif // MG_34
#ifndef AL_34
    #define AL_34
#endif // AL_34
#ifndef SI_34
    #define SI_34
#endif // SI_34
#ifndef P_34
    #define P_34
#endif // P_34
#ifndef S_34
    #define S_34
#endif // S_34
#ifndef CL_34
    #define CL_34
#endif // CL_34
#ifndef AR_34
    #define AR_34
#endif // AR_34
#ifndef K_34
    #define K_34
#endif // K_34
#ifndef CA_34
    #define CA_34
#endif // CA_34
#ifndef NA_35
    #define NA_35
#endif // NA_35
#ifndef MG_35
    #define MG_35
#endif // MG_35
#ifndef AL_35
    #define AL_35
#endif // AL_35
#ifndef SI_35
    #define SI_35
#endif // SI_35
#ifndef P_35
    #define P_35
#endif // P_35
#ifndef S_35
    #define S_35
#endif // S_35
#ifndef CL_35
    #define CL_35
#endif // CL_35
#ifndef AR_35
    #define AR_35
#endif // AR_35
#ifndef K_35
    #define K_35
#endif // K_35
#ifndef CA_35
    #define CA_35
#endif // CA_35
#ifndef SC_35
    #define SC_35
#endif // SC_35
#ifndef NA_36
    #define NA_36
#endif // NA_36
#ifndef MG_36
    #define MG_36
#endif // MG_36
#ifndef AL_36
    #define AL_36
#endif // AL_36
#ifndef SI_36
    #define SI_36
#endif // SI_36
#ifndef P_36
    #define P_36
#endif // P_36
#ifndef S_36
    #define S_36
#endif // S_36
#ifndef CL_36
    #define CL_36
#endif // CL_36
#ifndef AR_36
    #define AR_36
#endif // AR_36
#ifndef K_36
    #define K_36
#endif // K_36
#ifndef CA_36
    #define CA_36
#endif // CA_36
#ifndef SC_36
    #define SC_36
#endif // SC_36
#ifndef NA_37
    #define NA_37
#endif // NA_37
#ifndef MG_37
    #define MG_37
#endif // MG_37
#ifndef AL_37
    #define AL_37
#endif // AL_37
#ifndef SI_37
    #define SI_37
#endif // SI_37
#ifndef P_37
    #define P_37
#endif // P_37
#ifndef S_37
    #define S_37
#endif // S_37
#ifndef CL_37
    #define CL_37
#endif // CL_37
#ifndef AR_37
    #define AR_37
#endif // AR_37
#ifndef K_37
    #define K_37
#endif // K_37
#ifndef CA_37
    #define CA_37
#endif // CA_37
#ifndef SC_37
    #define SC_37
#endif // SC_37
#ifndef TI_37
    #define TI_37
#endif // TI_37
#ifndef NA_38
    #define NA_38
#endif // NA_38
#ifndef MG_38
    #define MG_38
#endif // MG_38
#ifndef AL_38
    #define AL_38
#endif // AL_38
#ifndef SI_38
    #define SI_38
#endif // SI_38
#ifndef P_38
    #define P_38
#endif // P_38
#ifndef S_38
    #define S_38
#endif // S_38
#ifndef CL_38
    #define CL_38
#endif // CL_38
#ifndef AR_38
    #define AR_38
#endif // AR_38
#ifndef K_38
    #define K_38
#endif // K_38
#ifndef CA_38
    #define CA_38
#endif // CA_38
#ifndef SC_38
    #define SC_38
#endif // SC_38
#ifndef TI_38
    #define TI_38
#endif // TI_38
#ifndef NA_39
    #define NA_39
#endif // NA_39
#ifndef MG_39
    #define MG_39
#endif // MG_39
#ifndef AL_39
    #define AL_39
#endif // AL_39
#ifndef SI_39
    #define SI_39
#endif // SI_39
#ifndef P_39
    #define P_39
#endif // P_39
#ifndef S_39
    #define S_39
#endif // S_39
#ifndef CL_39
    #define CL_39
#endif // CL_39
#ifndef AR_39
    #define AR_39
#endif // AR_39
#ifndef K_39
    #define K_39
#endif // K_39
#ifndef CA_39
    #define CA_39
#endif // CA_39
#ifndef SC_39
    #define SC_39
#endif // SC_39
#ifndef TI_39
    #define TI_39
#endif // TI_39
#ifndef V_39
    #define V_39
#endif // V_39
#ifndef MG_40
    #define MG_40
#endif // MG_40
#ifndef AL_40
    #define AL_40
#endif // AL_40
#ifndef SI_40
    #define SI_40
#endif // SI_40
#ifndef P_40
    #define P_40
#endif // P_40
#ifndef S_40
    #define S_40
#endif // S_40
#ifndef CL_40
    #define CL_40
#endif // CL_40
#ifndef AR_40
    #define AR_40
#endif // AR_40
#ifndef K_40
    #define K_40
#endif // K_40
#ifndef CA_40
    #define CA_40
#endif // CA_40
#ifndef SC_40
    #define SC_40
#endif // SC_40
#ifndef TI_40
    #define TI_40
#endif // TI_40
#ifndef V_40
    #define V_40
#endif // V_40
#ifndef MG_41
    #define MG_41
#endif // MG_41
#ifndef AL_41
    #define AL_41
#endif // AL_41
#ifndef SI_41
    #define SI_41
#endif // SI_41
#ifndef P_41
    #define P_41
#endif // P_41
#ifndef S_41
    #define S_41
#endif // S_41
#ifndef CL_41
    #define CL_41
#endif // CL_41
#ifndef AR_41
    #define AR_41
#endif // AR_41
#ifndef K_41
    #define K_41
#endif // K_41
#ifndef CA_41
    #define CA_41
#endif // CA_41
#ifndef SC_41
    #define SC_41
#endif // SC_41
#ifndef TI_41
    #define TI_41
#endif // TI_41
#ifndef V_41
    #define V_41
#endif // V_41
#ifndef CR_41
    #define CR_41
#endif // CR_41
#ifndef AL_42
    #define AL_42
#endif // AL_42
#ifndef SI_42
    #define SI_42
#endif // SI_42
#ifndef P_42
    #define P_42
#endif // P_42
#ifndef S_42
    #define S_42
#endif // S_42
#ifndef CL_42
    #define CL_42
#endif // CL_42
#ifndef AR_42
    #define AR_42
#endif // AR_42
#ifndef K_42
    #define K_42
#endif // K_42
#ifndef CA_42
    #define CA_42
#endif // CA_42
#ifndef SC_42
    #define SC_42
#endif // SC_42
#ifndef TI_42
    #define TI_42
#endif // TI_42
#ifndef V_42
    #define V_42
#endif // V_42
#ifndef CR_42
    #define CR_42
#endif // CR_42
#ifndef AL_43
    #define AL_43
#endif // AL_43
#ifndef SI_43
    #define SI_43
#endif // SI_43
#ifndef P_43
    #define P_43
#endif // P_43
#ifndef S_43
    #define S_43
#endif // S_43
#ifndef CL_43
    #define CL_43
#endif // CL_43
#ifndef AR_43
    #define AR_43
#endif // AR_43
#ifndef K_43
    #define K_43
#endif // K_43
#ifndef CA_43
    #define CA_43
#endif // CA_43
#ifndef SC_43
    #define SC_43
#endif // SC_43
#ifndef TI_43
    #define TI_43
#endif // TI_43
#ifndef V_43
    #define V_43
#endif // V_43
#ifndef CR_43
    #define CR_43
#endif // CR_43
#ifndef MN_43
    #define MN_43
#endif // MN_43
#ifndef SI_44
    #define SI_44
#endif // SI_44
#ifndef P_44
    #define P_44
#endif // P_44
#ifndef S_44
    #define S_44
#endif // S_44
#ifndef CL_44
    #define CL_44
#endif // CL_44
#ifndef AR_44
    #define AR_44
#endif // AR_44
#ifndef K_44
    #define K_44
#endif // K_44
#ifndef CA_44
    #define CA_44
#endif // CA_44
#ifndef SC_44
    #define SC_44
#endif // SC_44
#ifndef TI_44
    #define TI_44
#endif // TI_44
#ifndef V_44
    #define V_44
#endif // V_44
#ifndef CR_44
    #define CR_44
#endif // CR_44
#ifndef MN_44
    #define MN_44
#endif // MN_44
#ifndef SI_45
    #define SI_45
#endif // SI_45
#ifndef P_45
    #define P_45
#endif // P_45
#ifndef S_45
    #define S_45
#endif // S_45
#ifndef CL_45
    #define CL_45
#endif // CL_45
#ifndef AR_45
    #define AR_45
#endif // AR_45
#ifndef K_45
    #define K_45
#endif // K_45
#ifndef CA_45
    #define CA_45
#endif // CA_45
#ifndef SC_45
    #define SC_45
#endif // SC_45
#ifndef TI_45
    #define TI_45
#endif // TI_45
#ifndef V_45
    #define V_45
#endif // V_45
#ifndef CR_45
    #define CR_45
#endif // CR_45
#ifndef MN_45
    #define MN_45
#endif // MN_45
#ifndef FE_45
    #define FE_45
#endif // FE_45
#ifndef P_46
    #define P_46
#endif // P_46
#ifndef S_46
    #define S_46
#endif // S_46
#ifndef CL_46
    #define CL_46
#endif // CL_46
#ifndef AR_46
    #define AR_46
#endif // AR_46
#ifndef K_46
    #define K_46
#endif // K_46
#ifndef CA_46
    #define CA_46
#endif // CA_46
#ifndef SC_46
    #define SC_46
#endif // SC_46
#ifndef TI_46
    #define TI_46
#endif // TI_46
#ifndef V_46
    #define V_46
#endif // V_46
#ifndef CR_46
    #define CR_46
#endif // CR_46
#ifndef MN_46
    #define MN_46
#endif // MN_46
#ifndef FE_46
    #define FE_46
#endif // FE_46
#ifndef P_47
    #define P_47
#endif // P_47
#ifndef S_47
    #define S_47
#endif // S_47
#ifndef CL_47
    #define CL_47
#endif // CL_47
#ifndef AR_47
    #define AR_47
#endif // AR_47
#ifndef K_47
    #define K_47
#endif // K_47
#ifndef CA_47
    #define CA_47
#endif // CA_47
#ifndef SC_47
    #define SC_47
#endif // SC_47
#ifndef TI_47
    #define TI_47
#endif // TI_47
#ifndef V_47
    #define V_47
#endif // V_47
#ifndef CR_47
    #define CR_47
#endif // CR_47
#ifndef MN_47
    #define MN_47
#endif // MN_47
#ifndef FE_47
    #define FE_47
#endif // FE_47
#ifndef CO_47
    #define CO_47
#endif // CO_47
#ifndef S_48
    #define S_48
#endif // S_48
#ifndef CL_48
    #define CL_48
#endif // CL_48
#ifndef AR_48
    #define AR_48
#endif // AR_48
#ifndef K_48
    #define K_48
#endif // K_48
#ifndef CA_48
    #define CA_48
#endif // CA_48
#ifndef SC_48
    #define SC_48
#endif // SC_48
#ifndef TI_48
    #define TI_48
#endif // TI_48
#ifndef V_48
    #define V_48
#endif // V_48
#ifndef CR_48
    #define CR_48
#endif // CR_48
#ifndef MN_48
    #define MN_48
#endif // MN_48
#ifndef FE_48
    #define FE_48
#endif // FE_48
#ifndef CO_48
    #define CO_48
#endif // CO_48
#ifndef NI_48
    #define NI_48
#endif // NI_48
#ifndef S_49
    #define S_49
#endif // S_49
#ifndef CL_49
    #define CL_49
#endif // CL_49
#ifndef AR_49
    #define AR_49
#endif // AR_49
#ifndef K_49
    #define K_49
#endif // K_49
#ifndef CA_49
    #define CA_49
#endif // CA_49
#ifndef SC_49
    #define SC_49
#endif // SC_49
#ifndef TI_49
    #define TI_49
#endif // TI_49
#ifndef V_49
    #define V_49
#endif // V_49
#ifndef CR_49
    #define CR_49
#endif // CR_49
#ifndef MN_49
    #define MN_49
#endif // MN_49
#ifndef FE_49
    #define FE_49
#endif // FE_49
#ifndef CO_49
    #define CO_49
#endif // CO_49
#ifndef NI_49
    #define NI_49
#endif // NI_49
#ifndef CL_50
    #define CL_50
#endif // CL_50
#ifndef AR_50
    #define AR_50
#endif // AR_50
#ifndef K_50
    #define K_50
#endif // K_50
#ifndef CA_50
    #define CA_50
#endif // CA_50
#ifndef SC_50
    #define SC_50
#endif // SC_50
#ifndef TI_50
    #define TI_50
#endif // TI_50
#ifndef V_50
    #define V_50
#endif // V_50
#ifndef CR_50
    #define CR_50
#endif // CR_50
#ifndef MN_50
    #define MN_50
#endif // MN_50
#ifndef FE_50
    #define FE_50
#endif // FE_50
#ifndef CO_50
    #define CO_50
#endif // CO_50
#ifndef NI_50
    #define NI_50
#endif // NI_50
#ifndef CL_51
    #define CL_51
#endif // CL_51
#ifndef AR_51
    #define AR_51
#endif // AR_51
#ifndef K_51
    #define K_51
#endif // K_51
#ifndef CA_51
    #define CA_51
#endif // CA_51
#ifndef SC_51
    #define SC_51
#endif // SC_51
#ifndef TI_51
    #define TI_51
#endif // TI_51
#ifndef V_51
    #define V_51
#endif // V_51
#ifndef CR_51
    #define CR_51
#endif // CR_51
#ifndef MN_51
    #define MN_51
#endif // MN_51
#ifndef FE_51
    #define FE_51
#endif // FE_51
#ifndef CO_51
    #define CO_51
#endif // CO_51
#ifndef NI_51
    #define NI_51
#endif // NI_51
#ifndef CL_52
    #define CL_52
#endif // CL_52
#ifndef AR_52
    #define AR_52
#endif // AR_52
#ifndef K_52
    #define K_52
#endif // K_52
#ifndef CA_52
    #define CA_52
#endif // CA_52
#ifndef SC_52
    #define SC_52
#endif // SC_52
#ifndef TI_52
    #define TI_52
#endif // TI_52
#ifndef V_52
    #define V_52
#endif // V_52
#ifndef CR_52
    #define CR_52
#endif // CR_52
#ifndef MN_52
    #define MN_52
#endif // MN_52
#ifndef FE_52
    #define FE_52
#endif // FE_52
#ifndef CO_52
    #define CO_52
#endif // CO_52
#ifndef NI_52
    #define NI_52
#endif // NI_52
#ifndef CU_52
    #define CU_52
#endif // CU_52
#ifndef AR_53
    #define AR_53
#endif // AR_53
#ifndef K_53
    #define K_53
#endif // K_53
#ifndef CA_53
    #define CA_53
#endif // CA_53
#ifndef SC_53
    #define SC_53
#endif // SC_53
#ifndef TI_53
    #define TI_53
#endif // TI_53
#ifndef V_53
    #define V_53
#endif // V_53
#ifndef CR_53
    #define CR_53
#endif // CR_53
#ifndef MN_53
    #define MN_53
#endif // MN_53
#ifndef FE_53
    #define FE_53
#endif // FE_53
#ifndef CO_53
    #define CO_53
#endif // CO_53
#ifndef NI_53
    #define NI_53
#endif // NI_53
#ifndef CU_53
    #define CU_53
#endif // CU_53
#ifndef AR_54
    #define AR_54
#endif // AR_54
#ifndef K_54
    #define K_54
#endif // K_54
#ifndef CA_54
    #define CA_54
#endif // CA_54
#ifndef SC_54
    #define SC_54
#endif // SC_54
#ifndef TI_54
    #define TI_54
#endif // TI_54
#ifndef V_54
    #define V_54
#endif // V_54
#ifndef CR_54
    #define CR_54
#endif // CR_54
#ifndef MN_54
    #define MN_54
#endif // MN_54
#ifndef FE_54
    #define FE_54
#endif // FE_54
#ifndef CO_54
    #define CO_54
#endif // CO_54
#ifndef NI_54
    #define NI_54
#endif // NI_54
#ifndef CU_54
    #define CU_54
#endif // CU_54
#ifndef ZN_54
    #define ZN_54
#endif // ZN_54
#ifndef K_55
    #define K_55
#endif // K_55
#ifndef CA_55
    #define CA_55
#endif // CA_55
#ifndef SC_55
    #define SC_55
#endif // SC_55
#ifndef TI_55
    #define TI_55
#endif // TI_55
#ifndef V_55
    #define V_55
#endif // V_55
#ifndef CR_55
    #define CR_55
#endif // CR_55
#ifndef MN_55
    #define MN_55
#endif // MN_55
#ifndef FE_55
    #define FE_55
#endif // FE_55
#ifndef CO_55
    #define CO_55
#endif // CO_55
#ifndef NI_55
    #define NI_55
#endif // NI_55
#ifndef CU_55
    #define CU_55
#endif // CU_55
#ifndef ZN_55
    #define ZN_55
#endif // ZN_55
#ifndef K_56
    #define K_56
#endif // K_56
#ifndef CA_56
    #define CA_56
#endif // CA_56
#ifndef SC_56
    #define SC_56
#endif // SC_56
#ifndef TI_56
    #define TI_56
#endif // TI_56
#ifndef V_56
    #define V_56
#endif // V_56
#ifndef CR_56
    #define CR_56
#endif // CR_56
#ifndef MN_56
    #define MN_56
#endif // MN_56
#ifndef FE_56
    #define FE_56
#endif // FE_56
#ifndef CO_56
    #define CO_56
#endif // CO_56
#ifndef NI_56
    #define NI_56
#endif // NI_56
#ifndef CU_56
    #define CU_56
#endif // CU_56
#ifndef ZN_56
    #define ZN_56
#endif // ZN_56
#ifndef GA_56
    #define GA_56
#endif // GA_56
#ifndef K_57
    #define K_57
#endif // K_57
#ifndef CA_57
    #define CA_57
#endif // CA_57
#ifndef SC_57
    #define SC_57
#endif // SC_57
#ifndef TI_57
    #define TI_57
#endif // TI_57
#ifndef V_57
    #define V_57
#endif // V_57
#ifndef CR_57
    #define CR_57
#endif // CR_57
#ifndef MN_57
    #define MN_57
#endif // MN_57
#ifndef FE_57
    #define FE_57
#endif // FE_57
#ifndef CO_57
    #define CO_57
#endif // CO_57
#ifndef NI_57
    #define NI_57
#endif // NI_57
#ifndef CU_57
    #define CU_57
#endif // CU_57
#ifndef ZN_57
    #define ZN_57
#endif // ZN_57
#ifndef GA_57
    #define GA_57
#endif // GA_57
#ifndef K_58
    #define K_58
#endif // K_58
#ifndef CA_58
    #define CA_58
#endif // CA_58
#ifndef SC_58
    #define SC_58
#endif // SC_58
#ifndef TI_58
    #define TI_58
#endif // TI_58
#ifndef V_58
    #define V_58
#endif // V_58
#ifndef CR_58
    #define CR_58
#endif // CR_58
#ifndef MN_58
    #define MN_58
#endif // MN_58
#ifndef FE_58
    #define FE_58
#endif // FE_58
#ifndef CO_58
    #define CO_58
#endif // CO_58
#ifndef NI_58
    #define NI_58
#endif // NI_58
#ifndef CU_58
    #define CU_58
#endif // CU_58
#ifndef ZN_58
    #define ZN_58
#endif // ZN_58
#ifndef GA_58
    #define GA_58
#endif // GA_58
#ifndef GE_58
    #define GE_58
#endif // GE_58
#ifndef K_59
    #define K_59
#endif // K_59
#ifndef CA_59
    #define CA_59
#endif // CA_59
#ifndef SC_59
    #define SC_59
#endif // SC_59
#ifndef TI_59
    #define TI_59
#endif // TI_59
#ifndef V_59
    #define V_59
#endif // V_59
#ifndef CR_59
    #define CR_59
#endif // CR_59
#ifndef MN_59
    #define MN_59
#endif // MN_59
#ifndef FE_59
    #define FE_59
#endif // FE_59
#ifndef CO_59
    #define CO_59
#endif // CO_59
#ifndef NI_59
    #define NI_59
#endif // NI_59
#ifndef CU_59
    #define CU_59
#endif // CU_59
#ifndef ZN_59
    #define ZN_59
#endif // ZN_59
#ifndef GA_59
    #define GA_59
#endif // GA_59
#ifndef GE_59
    #define GE_59
#endif // GE_59
#ifndef CA_60
    #define CA_60
#endif // CA_60
#ifndef SC_60
    #define SC_60
#endif // SC_60
#ifndef TI_60
    #define TI_60
#endif // TI_60
#ifndef V_60
    #define V_60
#endif // V_60
#ifndef CR_60
    #define CR_60
#endif // CR_60
#ifndef MN_60
    #define MN_60
#endif // MN_60
#ifndef FE_60
    #define FE_60
#endif // FE_60
#ifndef CO_60
    #define CO_60
#endif // CO_60
#ifndef NI_60
    #define NI_60
#endif // NI_60
#ifndef CU_60
    #define CU_60
#endif // CU_60
#ifndef ZN_60
    #define ZN_60
#endif // ZN_60
#ifndef GA_60
    #define GA_60
#endif // GA_60
#ifndef GE_60
    #define GE_60
#endif // GE_60
#ifndef AS_60
    #define AS_60
#endif // AS_60
#ifndef CA_61
    #define CA_61
#endif // CA_61
#ifndef SC_61
    #define SC_61
#endif // SC_61
#ifndef TI_61
    #define TI_61
#endif // TI_61
#ifndef V_61
    #define V_61
#endif // V_61
#ifndef CR_61
    #define CR_61
#endif // CR_61
#ifndef MN_61
    #define MN_61
#endif // MN_61
#ifndef FE_61
    #define FE_61
#endif // FE_61
#ifndef CO_61
    #define CO_61
#endif // CO_61
#ifndef NI_61
    #define NI_61
#endif // NI_61
#ifndef CU_61
    #define CU_61
#endif // CU_61
#ifndef ZN_61
    #define ZN_61
#endif // ZN_61
#ifndef GA_61
    #define GA_61
#endif // GA_61
#ifndef GE_61
    #define GE_61
#endif // GE_61
#ifndef AS_61
    #define AS_61
#endif // AS_61
#ifndef SC_62
    #define SC_62
#endif // SC_62
#ifndef TI_62
    #define TI_62
#endif // TI_62
#ifndef V_62
    #define V_62
#endif // V_62
#ifndef CR_62
    #define CR_62
#endif // CR_62
#ifndef MN_62
    #define MN_62
#endif // MN_62
#ifndef FE_62
    #define FE_62
#endif // FE_62
#ifndef CO_62
    #define CO_62
#endif // CO_62
#ifndef NI_62
    #define NI_62
#endif // NI_62
#ifndef CU_62
    #define CU_62
#endif // CU_62
#ifndef ZN_62
    #define ZN_62
#endif // ZN_62
#ifndef GA_62
    #define GA_62
#endif // GA_62
#ifndef GE_62
    #define GE_62
#endif // GE_62
#ifndef AS_62
    #define AS_62
#endif // AS_62
#ifndef SC_63
    #define SC_63
#endif // SC_63
#ifndef TI_63
    #define TI_63
#endif // TI_63
#ifndef V_63
    #define V_63
#endif // V_63
#ifndef CR_63
    #define CR_63
#endif // CR_63
#ifndef MN_63
    #define MN_63
#endif // MN_63
#ifndef FE_63
    #define FE_63
#endif // FE_63
#ifndef CO_63
    #define CO_63
#endif // CO_63
#ifndef NI_63
    #define NI_63
#endif // NI_63
#ifndef CU_63
    #define CU_63
#endif // CU_63
#ifndef ZN_63
    #define ZN_63
#endif // ZN_63
#ifndef GA_63
    #define GA_63
#endif // GA_63
#ifndef GE_63
    #define GE_63
#endif // GE_63
#ifndef AS_63
    #define AS_63
#endif // AS_63
#ifndef SE_63
    #define SE_63
#endif // SE_63
#ifndef TI_64
    #define TI_64
#endif // TI_64
#ifndef V_64
    #define V_64
#endif // V_64
#ifndef CR_64
    #define CR_64
#endif // CR_64
#ifndef MN_64
    #define MN_64
#endif // MN_64
#ifndef FE_64
    #define FE_64
#endif // FE_64
#ifndef CO_64
    #define CO_64
#endif // CO_64
#ifndef NI_64
    #define NI_64
#endif // NI_64
#ifndef CU_64
    #define CU_64
#endif // CU_64
#ifndef ZN_64
    #define ZN_64
#endif // ZN_64
#ifndef GA_64
    #define GA_64
#endif // GA_64
#ifndef GE_64
    #define GE_64
#endif // GE_64
#ifndef AS_64
    #define AS_64
#endif // AS_64
#ifndef SE_64
    #define SE_64
#endif // SE_64
#ifndef TI_65
    #define TI_65
#endif // TI_65
#ifndef V_65
    #define V_65
#endif // V_65
#ifndef CR_65
    #define CR_65
#endif // CR_65
#ifndef MN_65
    #define MN_65
#endif // MN_65
#ifndef FE_65
    #define FE_65
#endif // FE_65
#ifndef CO_65
    #define CO_65
#endif // CO_65
#ifndef NI_65
    #define NI_65
#endif // NI_65
#ifndef CU_65
    #define CU_65
#endif // CU_65
#ifndef ZN_65
    #define ZN_65
#endif // ZN_65
#ifndef GA_65
    #define GA_65
#endif // GA_65
#ifndef GE_65
    #define GE_65
#endif // GE_65
#ifndef AS_65
    #define AS_65
#endif // AS_65
#ifndef SE_65
    #define SE_65
#endif // SE_65
#ifndef BR_65
    #define BR_65
#endif // BR_65
#ifndef V_66
    #define V_66
#endif // V_66
#ifndef CR_66
    #define CR_66
#endif // CR_66
#ifndef MN_66
    #define MN_66
#endif // MN_66
#ifndef FE_66
    #define FE_66
#endif // FE_66
#ifndef CO_66
    #define CO_66
#endif // CO_66
#ifndef NI_66
    #define NI_66
#endif // NI_66
#ifndef CU_66
    #define CU_66
#endif // CU_66
#ifndef ZN_66
    #define ZN_66
#endif // ZN_66
#ifndef GA_66
    #define GA_66
#endif // GA_66
#ifndef GE_66
    #define GE_66
#endif // GE_66
#ifndef AS_66
    #define AS_66
#endif // AS_66
#ifndef SE_66
    #define SE_66
#endif // SE_66
#ifndef BR_66
    #define BR_66
#endif // BR_66
#ifndef V_67
    #define V_67
#endif // V_67
#ifndef CR_67
    #define CR_67
#endif // CR_67
#ifndef MN_67
    #define MN_67
#endif // MN_67
#ifndef FE_67
    #define FE_67
#endif // FE_67
#ifndef CO_67
    #define CO_67
#endif // CO_67
#ifndef NI_67
    #define NI_67
#endif // NI_67
#ifndef CU_67
    #define CU_67
#endif // CU_67
#ifndef ZN_67
    #define ZN_67
#endif // ZN_67
#ifndef GA_67
    #define GA_67
#endif // GA_67
#ifndef GE_67
    #define GE_67
#endif // GE_67
#ifndef AS_67
    #define AS_67
#endif // AS_67
#ifndef SE_67
    #define SE_67
#endif // SE_67
#ifndef BR_67
    #define BR_67
#endif // BR_67
#ifndef KR_67
    #define KR_67
#endif // KR_67
#ifndef CR_68
    #define CR_68
#endif // CR_68
#ifndef MN_68
    #define MN_68
#endif // MN_68
#ifndef FE_68
    #define FE_68
#endif // FE_68
#ifndef CO_68
    #define CO_68
#endif // CO_68
#ifndef NI_68
    #define NI_68
#endif // NI_68
#ifndef CU_68
    #define CU_68
#endif // CU_68
#ifndef ZN_68
    #define ZN_68
#endif // ZN_68
#ifndef GA_68
    #define GA_68
#endif // GA_68
#ifndef GE_68
    #define GE_68
#endif // GE_68
#ifndef AS_68
    #define AS_68
#endif // AS_68
#ifndef SE_68
    #define SE_68
#endif // SE_68
#ifndef BR_68
    #define BR_68
#endif // BR_68
#ifndef KR_68
    #define KR_68
#endif // KR_68
#ifndef CR_69
    #define CR_69
#endif // CR_69
#ifndef MN_69
    #define MN_69
#endif // MN_69
#ifndef FE_69
    #define FE_69
#endif // FE_69
#ifndef CO_69
    #define CO_69
#endif // CO_69
#ifndef NI_69
    #define NI_69
#endif // NI_69
#ifndef CU_69
    #define CU_69
#endif // CU_69
#ifndef ZN_69
    #define ZN_69
#endif // ZN_69
#ifndef GA_69
    #define GA_69
#endif // GA_69
#ifndef GE_69
    #define GE_69
#endif // GE_69
#ifndef AS_69
    #define AS_69
#endif // AS_69
#ifndef SE_69
    #define SE_69
#endif // SE_69
#ifndef BR_69
    #define BR_69
#endif // BR_69
#ifndef KR_69
    #define KR_69
#endif // KR_69
#ifndef CR_70
    #define CR_70
#endif // CR_70
#ifndef MN_70
    #define MN_70
#endif // MN_70
#ifndef FE_70
    #define FE_70
#endif // FE_70
#ifndef CO_70
    #define CO_70
#endif // CO_70
#ifndef NI_70
    #define NI_70
#endif // NI_70
#ifndef CU_70
    #define CU_70
#endif // CU_70
#ifndef ZN_70
    #define ZN_70
#endif // ZN_70
#ifndef GA_70
    #define GA_70
#endif // GA_70
#ifndef GE_70
    #define GE_70
#endif // GE_70
#ifndef AS_70
    #define AS_70
#endif // AS_70
#ifndef SE_70
    #define SE_70
#endif // SE_70
#ifndef BR_70
    #define BR_70
#endif // BR_70
#ifndef KR_70
    #define KR_70
#endif // KR_70
#ifndef MN_71
    #define MN_71
#endif // MN_71
#ifndef FE_71
    #define FE_71
#endif // FE_71
#ifndef CO_71
    #define CO_71
#endif // CO_71
#ifndef NI_71
    #define NI_71
#endif // NI_71
#ifndef CU_71
    #define CU_71
#endif // CU_71
#ifndef ZN_71
    #define ZN_71
#endif // ZN_71
#ifndef GA_71
    #define GA_71
#endif // GA_71
#ifndef GE_71
    #define GE_71
#endif // GE_71
#ifndef AS_71
    #define AS_71
#endif // AS_71
#ifndef SE_71
    #define SE_71
#endif // SE_71
#ifndef BR_71
    #define BR_71
#endif // BR_71
#ifndef KR_71
    #define KR_71
#endif // KR_71
#ifndef RB_71
    #define RB_71
#endif // RB_71
#ifndef MN_72
    #define MN_72
#endif // MN_72
#ifndef FE_72
    #define FE_72
#endif // FE_72
#ifndef CO_72
    #define CO_72
#endif // CO_72
#ifndef NI_72
    #define NI_72
#endif // NI_72
#ifndef CU_72
    #define CU_72
#endif // CU_72
#ifndef ZN_72
    #define ZN_72
#endif // ZN_72
#ifndef GA_72
    #define GA_72
#endif // GA_72
#ifndef GE_72
    #define GE_72
#endif // GE_72
#ifndef AS_72
    #define AS_72
#endif // AS_72
#ifndef SE_72
    #define SE_72
#endif // SE_72
#ifndef BR_72
    #define BR_72
#endif // BR_72
#ifndef KR_72
    #define KR_72
#endif // KR_72
#ifndef RB_72
    #define RB_72
#endif // RB_72
#ifndef MN_73
    #define MN_73
#endif // MN_73
#ifndef FE_73
    #define FE_73
#endif // FE_73
#ifndef CO_73
    #define CO_73
#endif // CO_73
#ifndef NI_73
    #define NI_73
#endif // NI_73
#ifndef CU_73
    #define CU_73
#endif // CU_73
#ifndef ZN_73
    #define ZN_73
#endif // ZN_73
#ifndef GA_73
    #define GA_73
#endif // GA_73
#ifndef GE_73
    #define GE_73
#endif // GE_73
#ifndef AS_73
    #define AS_73
#endif // AS_73
#ifndef SE_73
    #define SE_73
#endif // SE_73
#ifndef BR_73
    #define BR_73
#endif // BR_73
#ifndef KR_73
    #define KR_73
#endif // KR_73
#ifndef RB_73
    #define RB_73
#endif // RB_73
#ifndef SR_73
    #define SR_73
#endif // SR_73
#ifndef FE_74
    #define FE_74
#endif // FE_74
#ifndef CO_74
    #define CO_74
#endif // CO_74
#ifndef NI_74
    #define NI_74
#endif // NI_74
#ifndef CU_74
    #define CU_74
#endif // CU_74
#ifndef ZN_74
    #define ZN_74
#endif // ZN_74
#ifndef GA_74
    #define GA_74
#endif // GA_74
#ifndef GE_74
    #define GE_74
#endif // GE_74
#ifndef AS_74
    #define AS_74
#endif // AS_74
#ifndef SE_74
    #define SE_74
#endif // SE_74
#ifndef BR_74
    #define BR_74
#endif // BR_74
#ifndef KR_74
    #define KR_74
#endif // KR_74
#ifndef RB_74
    #define RB_74
#endif // RB_74
#ifndef SR_74
    #define SR_74
#endif // SR_74
#ifndef FE_75
    #define FE_75
#endif // FE_75
#ifndef CO_75
    #define CO_75
#endif // CO_75
#ifndef NI_75
    #define NI_75
#endif // NI_75
#ifndef CU_75
    #define CU_75
#endif // CU_75
#ifndef ZN_75
    #define ZN_75
#endif // ZN_75
#ifndef GA_75
    #define GA_75
#endif // GA_75
#ifndef GE_75
    #define GE_75
#endif // GE_75
#ifndef AS_75
    #define AS_75
#endif // AS_75
#ifndef SE_75
    #define SE_75
#endif // SE_75
#ifndef BR_75
    #define BR_75
#endif // BR_75
#ifndef KR_75
    #define KR_75
#endif // KR_75
#ifndef RB_75
    #define RB_75
#endif // RB_75
#ifndef SR_75
    #define SR_75
#endif // SR_75
#ifndef Y_75
    #define Y_75
#endif // Y_75
#ifndef FE_76
    #define FE_76
#endif // FE_76
#ifndef CO_76
    #define CO_76
#endif // CO_76
#ifndef NI_76
    #define NI_76
#endif // NI_76
#ifndef CU_76
    #define CU_76
#endif // CU_76
#ifndef ZN_76
    #define ZN_76
#endif // ZN_76
#ifndef GA_76
    #define GA_76
#endif // GA_76
#ifndef GE_76
    #define GE_76
#endif // GE_76
#ifndef AS_76
    #define AS_76
#endif // AS_76
#ifndef SE_76
    #define SE_76
#endif // SE_76
#ifndef BR_76
    #define BR_76
#endif // BR_76
#ifndef KR_76
    #define KR_76
#endif // KR_76
#ifndef RB_76
    #define RB_76
#endif // RB_76
#ifndef SR_76
    #define SR_76
#endif // SR_76
#ifndef Y_76
    #define Y_76
#endif // Y_76
#ifndef CO_77
    #define CO_77
#endif // CO_77
#ifndef NI_77
    #define NI_77
#endif // NI_77
#ifndef CU_77
    #define CU_77
#endif // CU_77
#ifndef ZN_77
    #define ZN_77
#endif // ZN_77
#ifndef GA_77
    #define GA_77
#endif // GA_77
#ifndef GE_77
    #define GE_77
#endif // GE_77
#ifndef AS_77
    #define AS_77
#endif // AS_77
#ifndef SE_77
    #define SE_77
#endif // SE_77
#ifndef BR_77
    #define BR_77
#endif // BR_77
#ifndef KR_77
    #define KR_77
#endif // KR_77
#ifndef RB_77
    #define RB_77
#endif // RB_77
#ifndef SR_77
    #define SR_77
#endif // SR_77
#ifndef Y_77
    #define Y_77
#endif // Y_77
#ifndef ZR_77
    #define ZR_77
#endif // ZR_77
#ifndef CO_78
    #define CO_78
#endif // CO_78
#ifndef NI_78
    #define NI_78
#endif // NI_78
#ifndef CU_78
    #define CU_78
#endif // CU_78
#ifndef ZN_78
    #define ZN_78
#endif // ZN_78
#ifndef GA_78
    #define GA_78
#endif // GA_78
#ifndef GE_78
    #define GE_78
#endif // GE_78
#ifndef AS_78
    #define AS_78
#endif // AS_78
#ifndef SE_78
    #define SE_78
#endif // SE_78
#ifndef BR_78
    #define BR_78
#endif // BR_78
#ifndef KR_78
    #define KR_78
#endif // KR_78
#ifndef RB_78
    #define RB_78
#endif // RB_78
#ifndef SR_78
    #define SR_78
#endif // SR_78
#ifndef Y_78
    #define Y_78
#endif // Y_78
#ifndef ZR_78
    #define ZR_78
#endif // ZR_78
#ifndef NI_79
    #define NI_79
#endif // NI_79
#ifndef CU_79
    #define CU_79
#endif // CU_79
#ifndef ZN_79
    #define ZN_79
#endif // ZN_79
#ifndef GA_79
    #define GA_79
#endif // GA_79
#ifndef GE_79
    #define GE_79
#endif // GE_79
#ifndef AS_79
    #define AS_79
#endif // AS_79
#ifndef SE_79
    #define SE_79
#endif // SE_79
#ifndef BR_79
    #define BR_79
#endif // BR_79
#ifndef KR_79
    #define KR_79
#endif // KR_79
#ifndef RB_79
    #define RB_79
#endif // RB_79
#ifndef SR_79
    #define SR_79
#endif // SR_79
#ifndef Y_79
    #define Y_79
#endif // Y_79
#ifndef ZR_79
    #define ZR_79
#endif // ZR_79
#ifndef NB_79
    #define NB_79
#endif // NB_79
#ifndef NI_80
    #define NI_80
#endif // NI_80
#ifndef CU_80
    #define CU_80
#endif // CU_80
#ifndef ZN_80
    #define ZN_80
#endif // ZN_80
#ifndef GA_80
    #define GA_80
#endif // GA_80
#ifndef GE_80
    #define GE_80
#endif // GE_80
#ifndef AS_80
    #define AS_80
#endif // AS_80
#ifndef SE_80
    #define SE_80
#endif // SE_80
#ifndef BR_80
    #define BR_80
#endif // BR_80
#ifndef KR_80
    #define KR_80
#endif // KR_80
#ifndef RB_80
    #define RB_80
#endif // RB_80
#ifndef SR_80
    #define SR_80
#endif // SR_80
#ifndef Y_80
    #define Y_80
#endif // Y_80
#ifndef ZR_80
    #define ZR_80
#endif // ZR_80
#ifndef NB_80
    #define NB_80
#endif // NB_80
#ifndef NI_81
    #define NI_81
#endif // NI_81
#ifndef CU_81
    #define CU_81
#endif // CU_81
#ifndef ZN_81
    #define ZN_81
#endif // ZN_81
#ifndef GA_81
    #define GA_81
#endif // GA_81
#ifndef GE_81
    #define GE_81
#endif // GE_81
#ifndef AS_81
    #define AS_81
#endif // AS_81
#ifndef SE_81
    #define SE_81
#endif // SE_81
#ifndef BR_81
    #define BR_81
#endif // BR_81
#ifndef KR_81
    #define KR_81
#endif // KR_81
#ifndef RB_81
    #define RB_81
#endif // RB_81
#ifndef SR_81
    #define SR_81
#endif // SR_81
#ifndef Y_81
    #define Y_81
#endif // Y_81
#ifndef ZR_81
    #define ZR_81
#endif // ZR_81
#ifndef NB_81
    #define NB_81
#endif // NB_81
#ifndef MO_81
    #define MO_81
#endif // MO_81
#ifndef NI_82
    #define NI_82
#endif // NI_82
#ifndef CU_82
    #define CU_82
#endif // CU_82
#ifndef ZN_82
    #define ZN_82
#endif // ZN_82
#ifndef GA_82
    #define GA_82
#endif // GA_82
#ifndef GE_82
    #define GE_82
#endif // GE_82
#ifndef AS_82
    #define AS_82
#endif // AS_82
#ifndef SE_82
    #define SE_82
#endif // SE_82
#ifndef BR_82
    #define BR_82
#endif // BR_82
#ifndef KR_82
    #define KR_82
#endif // KR_82
#ifndef RB_82
    #define RB_82
#endif // RB_82
#ifndef SR_82
    #define SR_82
#endif // SR_82
#ifndef Y_82
    #define Y_82
#endif // Y_82
#ifndef ZR_82
    #define ZR_82
#endif // ZR_82
#ifndef NB_82
    #define NB_82
#endif // NB_82
#ifndef MO_82
    #define MO_82
#endif // MO_82
#ifndef CU_83
    #define CU_83
#endif // CU_83
#ifndef ZN_83
    #define ZN_83
#endif // ZN_83
#ifndef GA_83
    #define GA_83
#endif // GA_83
#ifndef GE_83
    #define GE_83
#endif // GE_83
#ifndef AS_83
    #define AS_83
#endif // AS_83
#ifndef SE_83
    #define SE_83
#endif // SE_83
#ifndef BR_83
    #define BR_83
#endif // BR_83
#ifndef KR_83
    #define KR_83
#endif // KR_83
#ifndef RB_83
    #define RB_83
#endif // RB_83
#ifndef SR_83
    #define SR_83
#endif // SR_83
#ifndef Y_83
    #define Y_83
#endif // Y_83
#ifndef ZR_83
    #define ZR_83
#endif // ZR_83
#ifndef NB_83
    #define NB_83
#endif // NB_83
#ifndef MO_83
    #define MO_83
#endif // MO_83
#ifndef TC_83
    #define TC_83
#endif // TC_83
#ifndef CU_84
    #define CU_84
#endif // CU_84
#ifndef ZN_84
    #define ZN_84
#endif // ZN_84
#ifndef GA_84
    #define GA_84
#endif // GA_84
#ifndef GE_84
    #define GE_84
#endif // GE_84
#ifndef AS_84
    #define AS_84
#endif // AS_84
#ifndef SE_84
    #define SE_84
#endif // SE_84
#ifndef BR_84
    #define BR_84
#endif // BR_84
#ifndef KR_84
    #define KR_84
#endif // KR_84
#ifndef RB_84
    #define RB_84
#endif // RB_84
#ifndef SR_84
    #define SR_84
#endif // SR_84
#ifndef Y_84
    #define Y_84
#endif // Y_84
#ifndef ZR_84
    #define ZR_84
#endif // ZR_84
#ifndef NB_84
    #define NB_84
#endif // NB_84
#ifndef MO_84
    #define MO_84
#endif // MO_84
#ifndef TC_84
    #define TC_84
#endif // TC_84
#ifndef ZN_85
    #define ZN_85
#endif // ZN_85
#ifndef GA_85
    #define GA_85
#endif // GA_85
#ifndef GE_85
    #define GE_85
#endif // GE_85
#ifndef AS_85
    #define AS_85
#endif // AS_85
#ifndef SE_85
    #define SE_85
#endif // SE_85
#ifndef BR_85
    #define BR_85
#endif // BR_85
#ifndef KR_85
    #define KR_85
#endif // KR_85
#ifndef RB_85
    #define RB_85
#endif // RB_85
#ifndef SR_85
    #define SR_85
#endif // SR_85
#ifndef Y_85
    #define Y_85
#endif // Y_85
#ifndef ZR_85
    #define ZR_85
#endif // ZR_85
#ifndef NB_85
    #define NB_85
#endif // NB_85
#ifndef MO_85
    #define MO_85
#endif // MO_85
#ifndef TC_85
    #define TC_85
#endif // TC_85
#ifndef RU_85
    #define RU_85
#endif // RU_85
#ifndef ZN_86
    #define ZN_86
#endif // ZN_86
#ifndef GA_86
    #define GA_86
#endif // GA_86
#ifndef GE_86
    #define GE_86
#endif // GE_86
#ifndef AS_86
    #define AS_86
#endif // AS_86
#ifndef SE_86
    #define SE_86
#endif // SE_86
#ifndef BR_86
    #define BR_86
#endif // BR_86
#ifndef KR_86
    #define KR_86
#endif // KR_86
#ifndef RB_86
    #define RB_86
#endif // RB_86
#ifndef SR_86
    #define SR_86
#endif // SR_86
#ifndef Y_86
    #define Y_86
#endif // Y_86
#ifndef ZR_86
    #define ZR_86
#endif // ZR_86
#ifndef NB_86
    #define NB_86
#endif // NB_86
#ifndef MO_86
    #define MO_86
#endif // MO_86
#ifndef TC_86
    #define TC_86
#endif // TC_86
#ifndef RU_86
    #define RU_86
#endif // RU_86
#ifndef GA_87
    #define GA_87
#endif // GA_87
#ifndef GE_87
    #define GE_87
#endif // GE_87
#ifndef AS_87
    #define AS_87
#endif // AS_87
#ifndef SE_87
    #define SE_87
#endif // SE_87
#ifndef BR_87
    #define BR_87
#endif // BR_87
#ifndef KR_87
    #define KR_87
#endif // KR_87
#ifndef RB_87
    #define RB_87
#endif // RB_87
#ifndef SR_87
    #define SR_87
#endif // SR_87
#ifndef Y_87
    #define Y_87
#endif // Y_87
#ifndef ZR_87
    #define ZR_87
#endif // ZR_87
#ifndef NB_87
    #define NB_87
#endif // NB_87
#ifndef MO_87
    #define MO_87
#endif // MO_87
#ifndef TC_87
    #define TC_87
#endif // TC_87
#ifndef RU_87
    #define RU_87
#endif // RU_87
#ifndef GA_88
    #define GA_88
#endif // GA_88
#ifndef GE_88
    #define GE_88
#endif // GE_88
#ifndef AS_88
    #define AS_88
#endif // AS_88
#ifndef SE_88
    #define SE_88
#endif // SE_88
#ifndef BR_88
    #define BR_88
#endif // BR_88
#ifndef KR_88
    #define KR_88
#endif // KR_88
#ifndef RB_88
    #define RB_88
#endif // RB_88
#ifndef SR_88
    #define SR_88
#endif // SR_88
#ifndef Y_88
    #define Y_88
#endif // Y_88
#ifndef ZR_88
    #define ZR_88
#endif // ZR_88
#ifndef NB_88
    #define NB_88
#endif // NB_88
#ifndef MO_88
    #define MO_88
#endif // MO_88
#ifndef TC_88
    #define TC_88
#endif // TC_88
#ifndef RU_88
    #define RU_88
#endif // RU_88
#ifndef RH_88
    #define RH_88
#endif // RH_88
#ifndef GE_89
    #define GE_89
#endif // GE_89
#ifndef AS_89
    #define AS_89
#endif // AS_89
#ifndef SE_89
    #define SE_89
#endif // SE_89
#ifndef BR_89
    #define BR_89
#endif // BR_89
#ifndef KR_89
    #define KR_89
#endif // KR_89
#ifndef RB_89
    #define RB_89
#endif // RB_89
#ifndef SR_89
    #define SR_89
#endif // SR_89
#ifndef Y_89
    #define Y_89
#endif // Y_89
#ifndef ZR_89
    #define ZR_89
#endif // ZR_89
#ifndef NB_89
    #define NB_89
#endif // NB_89
#ifndef MO_89
    #define MO_89
#endif // MO_89
#ifndef TC_89
    #define TC_89
#endif // TC_89
#ifndef RU_89
    #define RU_89
#endif // RU_89
#ifndef RH_89
    #define RH_89
#endif // RH_89
#ifndef GE_90
    #define GE_90
#endif // GE_90
#ifndef AS_90
    #define AS_90
#endif // AS_90
#ifndef SE_90
    #define SE_90
#endif // SE_90
#ifndef BR_90
    #define BR_90
#endif // BR_90
#ifndef KR_90
    #define KR_90
#endif // KR_90
#ifndef RB_90
    #define RB_90
#endif // RB_90
#ifndef SR_90
    #define SR_90
#endif // SR_90
#ifndef Y_90
    #define Y_90
#endif // Y_90
#ifndef ZR_90
    #define ZR_90
#endif // ZR_90
#ifndef NB_90
    #define NB_90
#endif // NB_90
#ifndef MO_90
    #define MO_90
#endif // MO_90
#ifndef TC_90
    #define TC_90
#endif // TC_90
#ifndef RU_90
    #define RU_90
#endif // RU_90
#ifndef RH_90
    #define RH_90
#endif // RH_90
#ifndef PD_90
    #define PD_90
#endif // PD_90
#ifndef AS_91
    #define AS_91
#endif // AS_91
#ifndef SE_91
    #define SE_91
#endif // SE_91
#ifndef BR_91
    #define BR_91
#endif // BR_91
#ifndef KR_91
    #define KR_91
#endif // KR_91
#ifndef RB_91
    #define RB_91
#endif // RB_91
#ifndef SR_91
    #define SR_91
#endif // SR_91
#ifndef Y_91
    #define Y_91
#endif // Y_91
#ifndef ZR_91
    #define ZR_91
#endif // ZR_91
#ifndef NB_91
    #define NB_91
#endif // NB_91
#ifndef MO_91
    #define MO_91
#endif // MO_91
#ifndef TC_91
    #define TC_91
#endif // TC_91
#ifndef RU_91
    #define RU_91
#endif // RU_91
#ifndef RH_91
    #define RH_91
#endif // RH_91
#ifndef PD_91
    #define PD_91
#endif // PD_91
#ifndef AS_92
    #define AS_92
#endif // AS_92
#ifndef SE_92
    #define SE_92
#endif // SE_92
#ifndef BR_92
    #define BR_92
#endif // BR_92
#ifndef KR_92
    #define KR_92
#endif // KR_92
#ifndef RB_92
    #define RB_92
#endif // RB_92
#ifndef SR_92
    #define SR_92
#endif // SR_92
#ifndef Y_92
    #define Y_92
#endif // Y_92
#ifndef ZR_92
    #define ZR_92
#endif // ZR_92
#ifndef NB_92
    #define NB_92
#endif // NB_92
#ifndef MO_92
    #define MO_92
#endif // MO_92
#ifndef TC_92
    #define TC_92
#endif // TC_92
#ifndef RU_92
    #define RU_92
#endif // RU_92
#ifndef RH_92
    #define RH_92
#endif // RH_92
#ifndef PD_92
    #define PD_92
#endif // PD_92
#ifndef AG_92
    #define AG_92
#endif // AG_92
#ifndef SE_93
    #define SE_93
#endif // SE_93
#ifndef BR_93
    #define BR_93
#endif // BR_93
#ifndef KR_93
    #define KR_93
#endif // KR_93
#ifndef RB_93
    #define RB_93
#endif // RB_93
#ifndef SR_93
    #define SR_93
#endif // SR_93
#ifndef Y_93
    #define Y_93
#endif // Y_93
#ifndef ZR_93
    #define ZR_93
#endif // ZR_93
#ifndef NB_93
    #define NB_93
#endif // NB_93
#ifndef MO_93
    #define MO_93
#endif // MO_93
#ifndef TC_93
    #define TC_93
#endif // TC_93
#ifndef RU_93
    #define RU_93
#endif // RU_93
#ifndef RH_93
    #define RH_93
#endif // RH_93
#ifndef PD_93
    #define PD_93
#endif // PD_93
#ifndef AG_93
    #define AG_93
#endif // AG_93
#ifndef SE_94
    #define SE_94
#endif // SE_94
#ifndef BR_94
    #define BR_94
#endif // BR_94
#ifndef KR_94
    #define KR_94
#endif // KR_94
#ifndef RB_94
    #define RB_94
#endif // RB_94
#ifndef SR_94
    #define SR_94
#endif // SR_94
#ifndef Y_94
    #define Y_94
#endif // Y_94
#ifndef ZR_94
    #define ZR_94
#endif // ZR_94
#ifndef NB_94
    #define NB_94
#endif // NB_94
#ifndef MO_94
    #define MO_94
#endif // MO_94
#ifndef TC_94
    #define TC_94
#endif // TC_94
#ifndef RU_94
    #define RU_94
#endif // RU_94
#ifndef RH_94
    #define RH_94
#endif // RH_94
#ifndef PD_94
    #define PD_94
#endif // PD_94
#ifndef AG_94
    #define AG_94
#endif // AG_94
#ifndef CD_94
    #define CD_94
#endif // CD_94
#ifndef SE_95
    #define SE_95
#endif // SE_95
#ifndef BR_95
    #define BR_95
#endif // BR_95
#ifndef KR_95
    #define KR_95
#endif // KR_95
#ifndef RB_95
    #define RB_95
#endif // RB_95
#ifndef SR_95
    #define SR_95
#endif // SR_95
#ifndef Y_95
    #define Y_95
#endif // Y_95
#ifndef ZR_95
    #define ZR_95
#endif // ZR_95
#ifndef NB_95
    #define NB_95
#endif // NB_95
#ifndef MO_95
    #define MO_95
#endif // MO_95
#ifndef TC_95
    #define TC_95
#endif // TC_95
#ifndef RU_95
    #define RU_95
#endif // RU_95
#ifndef RH_95
    #define RH_95
#endif // RH_95
#ifndef PD_95
    #define PD_95
#endif // PD_95
#ifndef AG_95
    #define AG_95
#endif // AG_95
#ifndef CD_95
    #define CD_95
#endif // CD_95
#ifndef BR_96
    #define BR_96
#endif // BR_96
#ifndef KR_96
    #define KR_96
#endif // KR_96
#ifndef RB_96
    #define RB_96
#endif // RB_96
#ifndef SR_96
    #define SR_96
#endif // SR_96
#ifndef Y_96
    #define Y_96
#endif // Y_96
#ifndef ZR_96
    #define ZR_96
#endif // ZR_96
#ifndef NB_96
    #define NB_96
#endif // NB_96
#ifndef MO_96
    #define MO_96
#endif // MO_96
#ifndef TC_96
    #define TC_96
#endif // TC_96
#ifndef RU_96
    #define RU_96
#endif // RU_96
#ifndef RH_96
    #define RH_96
#endif // RH_96
#ifndef PD_96
    #define PD_96
#endif // PD_96
#ifndef AG_96
    #define AG_96
#endif // AG_96
#ifndef CD_96
    #define CD_96
#endif // CD_96
#ifndef IN_96
    #define IN_96
#endif // IN_96
#ifndef BR_97
    #define BR_97
#endif // BR_97
#ifndef KR_97
    #define KR_97
#endif // KR_97
#ifndef RB_97
    #define RB_97
#endif // RB_97
#ifndef SR_97
    #define SR_97
#endif // SR_97
#ifndef Y_97
    #define Y_97
#endif // Y_97
#ifndef ZR_97
    #define ZR_97
#endif // ZR_97
#ifndef NB_97
    #define NB_97
#endif // NB_97
#ifndef MO_97
    #define MO_97
#endif // MO_97
#ifndef TC_97
    #define TC_97
#endif // TC_97
#ifndef RU_97
    #define RU_97
#endif // RU_97
#ifndef RH_97
    #define RH_97
#endif // RH_97
#ifndef PD_97
    #define PD_97
#endif // PD_97
#ifndef AG_97
    #define AG_97
#endif // AG_97
#ifndef CD_97
    #define CD_97
#endif // CD_97
#ifndef IN_97
    #define IN_97
#endif // IN_97
#ifndef BR_98
    #define BR_98
#endif // BR_98
#ifndef KR_98
    #define KR_98
#endif // KR_98
#ifndef RB_98
    #define RB_98
#endif // RB_98
#ifndef SR_98
    #define SR_98
#endif // SR_98
#ifndef Y_98
    #define Y_98
#endif // Y_98
#ifndef ZR_98
    #define ZR_98
#endif // ZR_98
#ifndef NB_98
    #define NB_98
#endif // NB_98
#ifndef MO_98
    #define MO_98
#endif // MO_98
#ifndef TC_98
    #define TC_98
#endif // TC_98
#ifndef RU_98
    #define RU_98
#endif // RU_98
#ifndef RH_98
    #define RH_98
#endif // RH_98
#ifndef PD_98
    #define PD_98
#endif // PD_98
#ifndef AG_98
    #define AG_98
#endif // AG_98
#ifndef CD_98
    #define CD_98
#endif // CD_98
#ifndef IN_98
    #define IN_98
#endif // IN_98
#ifndef KR_99
    #define KR_99
#endif // KR_99
#ifndef RB_99
    #define RB_99
#endif // RB_99
#ifndef SR_99
    #define SR_99
#endif // SR_99
#ifndef Y_99
    #define Y_99
#endif // Y_99
#ifndef ZR_99
    #define ZR_99
#endif // ZR_99
#ifndef NB_99
    #define NB_99
#endif // NB_99
#ifndef MO_99
    #define MO_99
#endif // MO_99
#ifndef TC_99
    #define TC_99
#endif // TC_99
#ifndef RU_99
    #define RU_99
#endif // RU_99
#ifndef RH_99
    #define RH_99
#endif // RH_99
#ifndef PD_99
    #define PD_99
#endif // PD_99
#ifndef AG_99
    #define AG_99
#endif // AG_99
#ifndef CD_99
    #define CD_99
#endif // CD_99
#ifndef IN_99
    #define IN_99
#endif // IN_99
#ifndef SN_99
    #define SN_99
#endif // SN_99
#ifndef KR_100
    #define KR_100
#endif // KR_100
#ifndef RB_100
    #define RB_100
#endif // RB_100
#ifndef SR_100
    #define SR_100
#endif // SR_100
#ifndef Y_100
    #define Y_100
#endif // Y_100
#ifndef ZR_100
    #define ZR_100
#endif // ZR_100
#ifndef NB_100
    #define NB_100
#endif // NB_100
#ifndef MO_100
    #define MO_100
#endif // MO_100
#ifndef TC_100
    #define TC_100
#endif // TC_100
#ifndef RU_100
    #define RU_100
#endif // RU_100
#ifndef RH_100
    #define RH_100
#endif // RH_100
#ifndef PD_100
    #define PD_100
#endif // PD_100
#ifndef AG_100
    #define AG_100
#endif // AG_100
#ifndef CD_100
    #define CD_100
#endif // CD_100
#ifndef IN_100
    #define IN_100
#endif // IN_100
#ifndef SN_100
    #define SN_100
#endif // SN_100
#ifndef KR_101
    #define KR_101
#endif // KR_101
#ifndef RB_101
    #define RB_101
#endif // RB_101
#ifndef SR_101
    #define SR_101
#endif // SR_101
#ifndef Y_101
    #define Y_101
#endif // Y_101
#ifndef ZR_101
    #define ZR_101
#endif // ZR_101
#ifndef NB_101
    #define NB_101
#endif // NB_101
#ifndef MO_101
    #define MO_101
#endif // MO_101
#ifndef TC_101
    #define TC_101
#endif // TC_101
#ifndef RU_101
    #define RU_101
#endif // RU_101
#ifndef RH_101
    #define RH_101
#endif // RH_101
#ifndef PD_101
    #define PD_101
#endif // PD_101
#ifndef AG_101
    #define AG_101
#endif // AG_101
#ifndef CD_101
    #define CD_101
#endif // CD_101
#ifndef IN_101
    #define IN_101
#endif // IN_101
#ifndef SN_101
    #define SN_101
#endif // SN_101
#ifndef RB_102
    #define RB_102
#endif // RB_102
#ifndef SR_102
    #define SR_102
#endif // SR_102
#ifndef Y_102
    #define Y_102
#endif // Y_102
#ifndef ZR_102
    #define ZR_102
#endif // ZR_102
#ifndef NB_102
    #define NB_102
#endif // NB_102
#ifndef MO_102
    #define MO_102
#endif // MO_102
#ifndef TC_102
    #define TC_102
#endif // TC_102
#ifndef RU_102
    #define RU_102
#endif // RU_102
#ifndef RH_102
    #define RH_102
#endif // RH_102
#ifndef PD_102
    #define PD_102
#endif // PD_102
#ifndef AG_102
    #define AG_102
#endif // AG_102
#ifndef CD_102
    #define CD_102
#endif // CD_102
#ifndef IN_102
    #define IN_102
#endif // IN_102
#ifndef SN_102
    #define SN_102
#endif // SN_102
#ifndef SB_102
    #define SB_102
#endif // SB_102
#ifndef RB_103
    #define RB_103
#endif // RB_103
#ifndef SR_103
    #define SR_103
#endif // SR_103
#ifndef Y_103
    #define Y_103
#endif // Y_103
#ifndef ZR_103
    #define ZR_103
#endif // ZR_103
#ifndef NB_103
    #define NB_103
#endif // NB_103
#ifndef MO_103
    #define MO_103
#endif // MO_103
#ifndef TC_103
    #define TC_103
#endif // TC_103
#ifndef RU_103
    #define RU_103
#endif // RU_103
#ifndef RH_103
    #define RH_103
#endif // RH_103
#ifndef PD_103
    #define PD_103
#endif // PD_103
#ifndef AG_103
    #define AG_103
#endif // AG_103
#ifndef CD_103
    #define CD_103
#endif // CD_103
#ifndef IN_103
    #define IN_103
#endif // IN_103
#ifndef SN_103
    #define SN_103
#endif // SN_103
#ifndef SB_103
    #define SB_103
#endif // SB_103
#ifndef RB_104
    #define RB_104
#endif // RB_104
#ifndef SR_104
    #define SR_104
#endif // SR_104
#ifndef Y_104
    #define Y_104
#endif // Y_104
#ifndef ZR_104
    #define ZR_104
#endif // ZR_104
#ifndef NB_104
    #define NB_104
#endif // NB_104
#ifndef MO_104
    #define MO_104
#endif // MO_104
#ifndef TC_104
    #define TC_104
#endif // TC_104
#ifndef RU_104
    #define RU_104
#endif // RU_104
#ifndef RH_104
    #define RH_104
#endif // RH_104
#ifndef PD_104
    #define PD_104
#endif // PD_104
#ifndef AG_104
    #define AG_104
#endif // AG_104
#ifndef CD_104
    #define CD_104
#endif // CD_104
#ifndef IN_104
    #define IN_104
#endif // IN_104
#ifndef SN_104
    #define SN_104
#endif // SN_104
#ifndef SB_104
    #define SB_104
#endif // SB_104
#ifndef TE_104
    #define TE_104
#endif // TE_104
#ifndef SR_105
    #define SR_105
#endif // SR_105
#ifndef Y_105
    #define Y_105
#endif // Y_105
#ifndef ZR_105
    #define ZR_105
#endif // ZR_105
#ifndef NB_105
    #define NB_105
#endif // NB_105
#ifndef MO_105
    #define MO_105
#endif // MO_105
#ifndef TC_105
    #define TC_105
#endif // TC_105
#ifndef RU_105
    #define RU_105
#endif // RU_105
#ifndef RH_105
    #define RH_105
#endif // RH_105
#ifndef PD_105
    #define PD_105
#endif // PD_105
#ifndef AG_105
    #define AG_105
#endif // AG_105
#ifndef CD_105
    #define CD_105
#endif // CD_105
#ifndef IN_105
    #define IN_105
#endif // IN_105
#ifndef SN_105
    #define SN_105
#endif // SN_105
#ifndef SB_105
    #define SB_105
#endif // SB_105
#ifndef TE_105
    #define TE_105
#endif // TE_105
#ifndef SR_106
    #define SR_106
#endif // SR_106
#ifndef Y_106
    #define Y_106
#endif // Y_106
#ifndef ZR_106
    #define ZR_106
#endif // ZR_106
#ifndef NB_106
    #define NB_106
#endif // NB_106
#ifndef MO_106
    #define MO_106
#endif // MO_106
#ifndef TC_106
    #define TC_106
#endif // TC_106
#ifndef RU_106
    #define RU_106
#endif // RU_106
#ifndef RH_106
    #define RH_106
#endif // RH_106
#ifndef PD_106
    #define PD_106
#endif // PD_106
#ifndef AG_106
    #define AG_106
#endif // AG_106
#ifndef CD_106
    #define CD_106
#endif // CD_106
#ifndef IN_106
    #define IN_106
#endif // IN_106
#ifndef SN_106
    #define SN_106
#endif // SN_106
#ifndef SB_106
    #define SB_106
#endif // SB_106
#ifndef TE_106
    #define TE_106
#endif // TE_106
#ifndef I_106
    #define I_106
#endif // I_106
#ifndef SR_107
    #define SR_107
#endif // SR_107
#ifndef Y_107
    #define Y_107
#endif // Y_107
#ifndef ZR_107
    #define ZR_107
#endif // ZR_107
#ifndef NB_107
    #define NB_107
#endif // NB_107
#ifndef MO_107
    #define MO_107
#endif // MO_107
#ifndef TC_107
    #define TC_107
#endif // TC_107
#ifndef RU_107
    #define RU_107
#endif // RU_107
#ifndef RH_107
    #define RH_107
#endif // RH_107
#ifndef PD_107
    #define PD_107
#endif // PD_107
#ifndef AG_107
    #define AG_107
#endif // AG_107
#ifndef CD_107
    #define CD_107
#endif // CD_107
#ifndef IN_107
    #define IN_107
#endif // IN_107
#ifndef SN_107
    #define SN_107
#endif // SN_107
#ifndef SB_107
    #define SB_107
#endif // SB_107
#ifndef TE_107
    #define TE_107
#endif // TE_107
#ifndef I_107
    #define I_107
#endif // I_107
#ifndef Y_108
    #define Y_108
#endif // Y_108
#ifndef ZR_108
    #define ZR_108
#endif // ZR_108
#ifndef NB_108
    #define NB_108
#endif // NB_108
#ifndef MO_108
    #define MO_108
#endif // MO_108
#ifndef TC_108
    #define TC_108
#endif // TC_108
#ifndef RU_108
    #define RU_108
#endif // RU_108
#ifndef RH_108
    #define RH_108
#endif // RH_108
#ifndef PD_108
    #define PD_108
#endif // PD_108
#ifndef AG_108
    #define AG_108
#endif // AG_108
#ifndef CD_108
    #define CD_108
#endif // CD_108
#ifndef IN_108
    #define IN_108
#endif // IN_108
#ifndef SN_108
    #define SN_108
#endif // SN_108
#ifndef SB_108
    #define SB_108
#endif // SB_108
#ifndef TE_108
    #define TE_108
#endif // TE_108
#ifndef I_108
    #define I_108
#endif // I_108
#ifndef XE_108
    #define XE_108
#endif // XE_108
#ifndef Y_109
    #define Y_109
#endif // Y_109
#ifndef ZR_109
    #define ZR_109
#endif // ZR_109
#ifndef NB_109
    #define NB_109
#endif // NB_109
#ifndef MO_109
    #define MO_109
#endif // MO_109
#ifndef TC_109
    #define TC_109
#endif // TC_109
#ifndef RU_109
    #define RU_109
#endif // RU_109
#ifndef RH_109
    #define RH_109
#endif // RH_109
#ifndef PD_109
    #define PD_109
#endif // PD_109
#ifndef AG_109
    #define AG_109
#endif // AG_109
#ifndef CD_109
    #define CD_109
#endif // CD_109
#ifndef IN_109
    #define IN_109
#endif // IN_109
#ifndef SN_109
    #define SN_109
#endif // SN_109
#ifndef SB_109
    #define SB_109
#endif // SB_109
#ifndef TE_109
    #define TE_109
#endif // TE_109
#ifndef I_109
    #define I_109
#endif // I_109
#ifndef XE_109
    #define XE_109
#endif // XE_109
#ifndef ZR_110
    #define ZR_110
#endif // ZR_110
#ifndef NB_110
    #define NB_110
#endif // NB_110
#ifndef MO_110
    #define MO_110
#endif // MO_110
#ifndef TC_110
    #define TC_110
#endif // TC_110
#ifndef RU_110
    #define RU_110
#endif // RU_110
#ifndef RH_110
    #define RH_110
#endif // RH_110
#ifndef PD_110
    #define PD_110
#endif // PD_110
#ifndef AG_110
    #define AG_110
#endif // AG_110
#ifndef CD_110
    #define CD_110
#endif // CD_110
#ifndef IN_110
    #define IN_110
#endif // IN_110
#ifndef SN_110
    #define SN_110
#endif // SN_110
#ifndef SB_110
    #define SB_110
#endif // SB_110
#ifndef TE_110
    #define TE_110
#endif // TE_110
#ifndef I_110
    #define I_110
#endif // I_110
#ifndef XE_110
    #define XE_110
#endif // XE_110
#ifndef ZR_111
    #define ZR_111
#endif // ZR_111
#ifndef NB_111
    #define NB_111
#endif // NB_111
#ifndef MO_111
    #define MO_111
#endif // MO_111
#ifndef TC_111
    #define TC_111
#endif // TC_111
#ifndef RU_111
    #define RU_111
#endif // RU_111
#ifndef RH_111
    #define RH_111
#endif // RH_111
#ifndef PD_111
    #define PD_111
#endif // PD_111
#ifndef AG_111
    #define AG_111
#endif // AG_111
#ifndef CD_111
    #define CD_111
#endif // CD_111
#ifndef IN_111
    #define IN_111
#endif // IN_111
#ifndef SN_111
    #define SN_111
#endif // SN_111
#ifndef SB_111
    #define SB_111
#endif // SB_111
#ifndef TE_111
    #define TE_111
#endif // TE_111
#ifndef I_111
    #define I_111
#endif // I_111
#ifndef XE_111
    #define XE_111
#endif // XE_111
#ifndef CS_111
    #define CS_111
#endif // CS_111
#ifndef ZR_112
    #define ZR_112
#endif // ZR_112
#ifndef NB_112
    #define NB_112
#endif // NB_112
#ifndef MO_112
    #define MO_112
#endif // MO_112
#ifndef TC_112
    #define TC_112
#endif // TC_112
#ifndef RU_112
    #define RU_112
#endif // RU_112
#ifndef RH_112
    #define RH_112
#endif // RH_112
#ifndef PD_112
    #define PD_112
#endif // PD_112
#ifndef AG_112
    #define AG_112
#endif // AG_112
#ifndef CD_112
    #define CD_112
#endif // CD_112
#ifndef IN_112
    #define IN_112
#endif // IN_112
#ifndef SN_112
    #define SN_112
#endif // SN_112
#ifndef SB_112
    #define SB_112
#endif // SB_112
#ifndef TE_112
    #define TE_112
#endif // TE_112
#ifndef I_112
    #define I_112
#endif // I_112
#ifndef XE_112
    #define XE_112
#endif // XE_112
#ifndef CS_112
    #define CS_112
#endif // CS_112
#ifndef ZR_113
    #define ZR_113
#endif // ZR_113
#ifndef NB_113
    #define NB_113
#endif // NB_113
#ifndef MO_113
    #define MO_113
#endif // MO_113
#ifndef TC_113
    #define TC_113
#endif // TC_113
#ifndef RU_113
    #define RU_113
#endif // RU_113
#ifndef RH_113
    #define RH_113
#endif // RH_113
#ifndef PD_113
    #define PD_113
#endif // PD_113
#ifndef AG_113
    #define AG_113
#endif // AG_113
#ifndef CD_113
    #define CD_113
#endif // CD_113
#ifndef IN_113
    #define IN_113
#endif // IN_113
#ifndef SN_113
    #define SN_113
#endif // SN_113
#ifndef SB_113
    #define SB_113
#endif // SB_113
#ifndef TE_113
    #define TE_113
#endif // TE_113
#ifndef I_113
    #define I_113
#endif // I_113
#ifndef XE_113
    #define XE_113
#endif // XE_113
#ifndef CS_113
    #define CS_113
#endif // CS_113
#ifndef BA_113
    #define BA_113
#endif // BA_113
#ifndef NB_114
    #define NB_114
#endif // NB_114
#ifndef MO_114
    #define MO_114
#endif // MO_114
#ifndef TC_114
    #define TC_114
#endif // TC_114
#ifndef RU_114
    #define RU_114
#endif // RU_114
#ifndef RH_114
    #define RH_114
#endif // RH_114
#ifndef PD_114
    #define PD_114
#endif // PD_114
#ifndef AG_114
    #define AG_114
#endif // AG_114
#ifndef CD_114
    #define CD_114
#endif // CD_114
#ifndef IN_114
    #define IN_114
#endif // IN_114
#ifndef SN_114
    #define SN_114
#endif // SN_114
#ifndef SB_114
    #define SB_114
#endif // SB_114
#ifndef TE_114
    #define TE_114
#endif // TE_114
#ifndef I_114
    #define I_114
#endif // I_114
#ifndef XE_114
    #define XE_114
#endif // XE_114
#ifndef CS_114
    #define CS_114
#endif // CS_114
#ifndef BA_114
    #define BA_114
#endif // BA_114
#ifndef NB_115
    #define NB_115
#endif // NB_115
#ifndef MO_115
    #define MO_115
#endif // MO_115
#ifndef TC_115
    #define TC_115
#endif // TC_115
#ifndef RU_115
    #define RU_115
#endif // RU_115
#ifndef RH_115
    #define RH_115
#endif // RH_115
#ifndef PD_115
    #define PD_115
#endif // PD_115
#ifndef AG_115
    #define AG_115
#endif // AG_115
#ifndef CD_115
    #define CD_115
#endif // CD_115
#ifndef IN_115
    #define IN_115
#endif // IN_115
#ifndef SN_115
    #define SN_115
#endif // SN_115
#ifndef SB_115
    #define SB_115
#endif // SB_115
#ifndef TE_115
    #define TE_115
#endif // TE_115
#ifndef I_115
    #define I_115
#endif // I_115
#ifndef XE_115
    #define XE_115
#endif // XE_115
#ifndef CS_115
    #define CS_115
#endif // CS_115
#ifndef BA_115
    #define BA_115
#endif // BA_115
#ifndef NB_116
    #define NB_116
#endif // NB_116
#ifndef MO_116
    #define MO_116
#endif // MO_116
#ifndef TC_116
    #define TC_116
#endif // TC_116
#ifndef RU_116
    #define RU_116
#endif // RU_116
#ifndef RH_116
    #define RH_116
#endif // RH_116
#ifndef PD_116
    #define PD_116
#endif // PD_116
#ifndef AG_116
    #define AG_116
#endif // AG_116
#ifndef CD_116
    #define CD_116
#endif // CD_116
#ifndef IN_116
    #define IN_116
#endif // IN_116
#ifndef SN_116
    #define SN_116
#endif // SN_116
#ifndef SB_116
    #define SB_116
#endif // SB_116
#ifndef TE_116
    #define TE_116
#endif // TE_116
#ifndef I_116
    #define I_116
#endif // I_116
#ifndef XE_116
    #define XE_116
#endif // XE_116
#ifndef CS_116
    #define CS_116
#endif // CS_116
#ifndef BA_116
    #define BA_116
#endif // BA_116
#ifndef LA_116
    #define LA_116
#endif // LA_116
#ifndef MO_117
    #define MO_117
#endif // MO_117
#ifndef TC_117
    #define TC_117
#endif // TC_117
#ifndef RU_117
    #define RU_117
#endif // RU_117
#ifndef RH_117
    #define RH_117
#endif // RH_117
#ifndef PD_117
    #define PD_117
#endif // PD_117
#ifndef AG_117
    #define AG_117
#endif // AG_117
#ifndef CD_117
    #define CD_117
#endif // CD_117
#ifndef IN_117
    #define IN_117
#endif // IN_117
#ifndef SN_117
    #define SN_117
#endif // SN_117
#ifndef SB_117
    #define SB_117
#endif // SB_117
#ifndef TE_117
    #define TE_117
#endif // TE_117
#ifndef I_117
    #define I_117
#endif // I_117
#ifndef XE_117
    #define XE_117
#endif // XE_117
#ifndef CS_117
    #define CS_117
#endif // CS_117
#ifndef BA_117
    #define BA_117
#endif // BA_117
#ifndef LA_117
    #define LA_117
#endif // LA_117
#ifndef MO_118
    #define MO_118
#endif // MO_118
#ifndef TC_118
    #define TC_118
#endif // TC_118
#ifndef RU_118
    #define RU_118
#endif // RU_118
#ifndef RH_118
    #define RH_118
#endif // RH_118
#ifndef PD_118
    #define PD_118
#endif // PD_118
#ifndef AG_118
    #define AG_118
#endif // AG_118
#ifndef CD_118
    #define CD_118
#endif // CD_118
#ifndef IN_118
    #define IN_118
#endif // IN_118
#ifndef SN_118
    #define SN_118
#endif // SN_118
#ifndef SB_118
    #define SB_118
#endif // SB_118
#ifndef TE_118
    #define TE_118
#endif // TE_118
#ifndef I_118
    #define I_118
#endif // I_118
#ifndef XE_118
    #define XE_118
#endif // XE_118
#ifndef CS_118
    #define CS_118
#endif // CS_118
#ifndef BA_118
    #define BA_118
#endif // BA_118
#ifndef LA_118
    #define LA_118
#endif // LA_118
#ifndef MO_119
    #define MO_119
#endif // MO_119
#ifndef TC_119
    #define TC_119
#endif // TC_119
#ifndef RU_119
    #define RU_119
#endif // RU_119
#ifndef RH_119
    #define RH_119
#endif // RH_119
#ifndef PD_119
    #define PD_119
#endif // PD_119
#ifndef AG_119
    #define AG_119
#endif // AG_119
#ifndef CD_119
    #define CD_119
#endif // CD_119
#ifndef IN_119
    #define IN_119
#endif // IN_119
#ifndef SN_119
    #define SN_119
#endif // SN_119
#ifndef SB_119
    #define SB_119
#endif // SB_119
#ifndef TE_119
    #define TE_119
#endif // TE_119
#ifndef I_119
    #define I_119
#endif // I_119
#ifndef XE_119
    #define XE_119
#endif // XE_119
#ifndef CS_119
    #define CS_119
#endif // CS_119
#ifndef BA_119
    #define BA_119
#endif // BA_119
#ifndef LA_119
    #define LA_119
#endif // LA_119
#ifndef CE_119
    #define CE_119
#endif // CE_119
#ifndef TC_120
    #define TC_120
#endif // TC_120
#ifndef RU_120
    #define RU_120
#endif // RU_120
#ifndef RH_120
    #define RH_120
#endif // RH_120
#ifndef PD_120
    #define PD_120
#endif // PD_120
#ifndef AG_120
    #define AG_120
#endif // AG_120
#ifndef CD_120
    #define CD_120
#endif // CD_120
#ifndef IN_120
    #define IN_120
#endif // IN_120
#ifndef SN_120
    #define SN_120
#endif // SN_120
#ifndef SB_120
    #define SB_120
#endif // SB_120
#ifndef TE_120
    #define TE_120
#endif // TE_120
#ifndef I_120
    #define I_120
#endif // I_120
#ifndef XE_120
    #define XE_120
#endif // XE_120
#ifndef CS_120
    #define CS_120
#endif // CS_120
#ifndef BA_120
    #define BA_120
#endif // BA_120
#ifndef LA_120
    #define LA_120
#endif // LA_120
#ifndef CE_120
    #define CE_120
#endif // CE_120
#ifndef TC_121
    #define TC_121
#endif // TC_121
#ifndef RU_121
    #define RU_121
#endif // RU_121
#ifndef RH_121
    #define RH_121
#endif // RH_121
#ifndef PD_121
    #define PD_121
#endif // PD_121
#ifndef AG_121
    #define AG_121
#endif // AG_121
#ifndef CD_121
    #define CD_121
#endif // CD_121
#ifndef IN_121
    #define IN_121
#endif // IN_121
#ifndef SN_121
    #define SN_121
#endif // SN_121
#ifndef SB_121
    #define SB_121
#endif // SB_121
#ifndef TE_121
    #define TE_121
#endif // TE_121
#ifndef I_121
    #define I_121
#endif // I_121
#ifndef XE_121
    #define XE_121
#endif // XE_121
#ifndef CS_121
    #define CS_121
#endif // CS_121
#ifndef BA_121
    #define BA_121
#endif // BA_121
#ifndef LA_121
    #define LA_121
#endif // LA_121
#ifndef CE_121
    #define CE_121
#endif // CE_121
#ifndef PR_121
    #define PR_121
#endif // PR_121
#ifndef TC_122
    #define TC_122
#endif // TC_122
#ifndef RU_122
    #define RU_122
#endif // RU_122
#ifndef RH_122
    #define RH_122
#endif // RH_122
#ifndef PD_122
    #define PD_122
#endif // PD_122
#ifndef AG_122
    #define AG_122
#endif // AG_122
#ifndef CD_122
    #define CD_122
#endif // CD_122
#ifndef IN_122
    #define IN_122
#endif // IN_122
#ifndef SN_122
    #define SN_122
#endif // SN_122
#ifndef SB_122
    #define SB_122
#endif // SB_122
#ifndef TE_122
    #define TE_122
#endif // TE_122
#ifndef I_122
    #define I_122
#endif // I_122
#ifndef XE_122
    #define XE_122
#endif // XE_122
#ifndef CS_122
    #define CS_122
#endif // CS_122
#ifndef BA_122
    #define BA_122
#endif // BA_122
#ifndef LA_122
    #define LA_122
#endif // LA_122
#ifndef CE_122
    #define CE_122
#endif // CE_122
#ifndef PR_122
    #define PR_122
#endif // PR_122
#ifndef RU_123
    #define RU_123
#endif // RU_123
#ifndef RH_123
    #define RH_123
#endif // RH_123
#ifndef PD_123
    #define PD_123
#endif // PD_123
#ifndef AG_123
    #define AG_123
#endif // AG_123
#ifndef CD_123
    #define CD_123
#endif // CD_123
#ifndef IN_123
    #define IN_123
#endif // IN_123
#ifndef SN_123
    #define SN_123
#endif // SN_123
#ifndef SB_123
    #define SB_123
#endif // SB_123
#ifndef TE_123
    #define TE_123
#endif // TE_123
#ifndef I_123
    #define I_123
#endif // I_123
#ifndef XE_123
    #define XE_123
#endif // XE_123
#ifndef CS_123
    #define CS_123
#endif // CS_123
#ifndef BA_123
    #define BA_123
#endif // BA_123
#ifndef LA_123
    #define LA_123
#endif // LA_123
#ifndef CE_123
    #define CE_123
#endif // CE_123
#ifndef PR_123
    #define PR_123
#endif // PR_123
#ifndef RU_124
    #define RU_124
#endif // RU_124
#ifndef RH_124
    #define RH_124
#endif // RH_124
#ifndef PD_124
    #define PD_124
#endif // PD_124
#ifndef AG_124
    #define AG_124
#endif // AG_124
#ifndef CD_124
    #define CD_124
#endif // CD_124
#ifndef IN_124
    #define IN_124
#endif // IN_124
#ifndef SN_124
    #define SN_124
#endif // SN_124
#ifndef SB_124
    #define SB_124
#endif // SB_124
#ifndef TE_124
    #define TE_124
#endif // TE_124
#ifndef I_124
    #define I_124
#endif // I_124
#ifndef XE_124
    #define XE_124
#endif // XE_124
#ifndef CS_124
    #define CS_124
#endif // CS_124
#ifndef BA_124
    #define BA_124
#endif // BA_124
#ifndef LA_124
    #define LA_124
#endif // LA_124
#ifndef CE_124
    #define CE_124
#endif // CE_124
#ifndef PR_124
    #define PR_124
#endif // PR_124
#ifndef ND_124
    #define ND_124
#endif // ND_124
#ifndef RU_125
    #define RU_125
#endif // RU_125
#ifndef RH_125
    #define RH_125
#endif // RH_125
#ifndef PD_125
    #define PD_125
#endif // PD_125
#ifndef AG_125
    #define AG_125
#endif // AG_125
#ifndef CD_125
    #define CD_125
#endif // CD_125
#ifndef IN_125
    #define IN_125
#endif // IN_125
#ifndef SN_125
    #define SN_125
#endif // SN_125
#ifndef SB_125
    #define SB_125
#endif // SB_125
#ifndef TE_125
    #define TE_125
#endif // TE_125
#ifndef I_125
    #define I_125
#endif // I_125
#ifndef XE_125
    #define XE_125
#endif // XE_125
#ifndef CS_125
    #define CS_125
#endif // CS_125
#ifndef BA_125
    #define BA_125
#endif // BA_125
#ifndef LA_125
    #define LA_125
#endif // LA_125
#ifndef CE_125
    #define CE_125
#endif // CE_125
#ifndef PR_125
    #define PR_125
#endif // PR_125
#ifndef ND_125
    #define ND_125
#endif // ND_125
#ifndef RH_126
    #define RH_126
#endif // RH_126
#ifndef PD_126
    #define PD_126
#endif // PD_126
#ifndef AG_126
    #define AG_126
#endif // AG_126
#ifndef CD_126
    #define CD_126
#endif // CD_126
#ifndef IN_126
    #define IN_126
#endif // IN_126
#ifndef SN_126
    #define SN_126
#endif // SN_126
#ifndef SB_126
    #define SB_126
#endif // SB_126
#ifndef TE_126
    #define TE_126
#endif // TE_126
#ifndef I_126
    #define I_126
#endif // I_126
#ifndef XE_126
    #define XE_126
#endif // XE_126
#ifndef CS_126
    #define CS_126
#endif // CS_126
#ifndef BA_126
    #define BA_126
#endif // BA_126
#ifndef LA_126
    #define LA_126
#endif // LA_126
#ifndef CE_126
    #define CE_126
#endif // CE_126
#ifndef PR_126
    #define PR_126
#endif // PR_126
#ifndef ND_126
    #define ND_126
#endif // ND_126
#ifndef PM_126
    #define PM_126
#endif // PM_126
#ifndef RH_127
    #define RH_127
#endif // RH_127
#ifndef PD_127
    #define PD_127
#endif // PD_127
#ifndef AG_127
    #define AG_127
#endif // AG_127
#ifndef CD_127
    #define CD_127
#endif // CD_127
#ifndef IN_127
    #define IN_127
#endif // IN_127
#ifndef SN_127
    #define SN_127
#endif // SN_127
#ifndef SB_127
    #define SB_127
#endif // SB_127
#ifndef TE_127
    #define TE_127
#endif // TE_127
#ifndef I_127
    #define I_127
#endif // I_127
#ifndef XE_127
    #define XE_127
#endif // XE_127
#ifndef CS_127
    #define CS_127
#endif // CS_127
#ifndef BA_127
    #define BA_127
#endif // BA_127
#ifndef LA_127
    #define LA_127
#endif // LA_127
#ifndef CE_127
    #define CE_127
#endif // CE_127
#ifndef PR_127
    #define PR_127
#endif // PR_127
#ifndef ND_127
    #define ND_127
#endif // ND_127
#ifndef PM_127
    #define PM_127
#endif // PM_127
#ifndef RH_128
    #define RH_128
#endif // RH_128
#ifndef PD_128
    #define PD_128
#endif // PD_128
#ifndef AG_128
    #define AG_128
#endif // AG_128
#ifndef CD_128
    #define CD_128
#endif // CD_128
#ifndef IN_128
    #define IN_128
#endif // IN_128
#ifndef SN_128
    #define SN_128
#endif // SN_128
#ifndef SB_128
    #define SB_128
#endif // SB_128
#ifndef TE_128
    #define TE_128
#endif // TE_128
#ifndef I_128
    #define I_128
#endif // I_128
#ifndef XE_128
    #define XE_128
#endif // XE_128
#ifndef CS_128
    #define CS_128
#endif // CS_128
#ifndef BA_128
    #define BA_128
#endif // BA_128
#ifndef LA_128
    #define LA_128
#endif // LA_128
#ifndef CE_128
    #define CE_128
#endif // CE_128
#ifndef PR_128
    #define PR_128
#endif // PR_128
#ifndef ND_128
    #define ND_128
#endif // ND_128
#ifndef PM_128
    #define PM_128
#endif // PM_128
#ifndef SM_128
    #define SM_128
#endif // SM_128
#ifndef PD_129
    #define PD_129
#endif // PD_129
#ifndef AG_129
    #define AG_129
#endif // AG_129
#ifndef CD_129
    #define CD_129
#endif // CD_129
#ifndef IN_129
    #define IN_129
#endif // IN_129
#ifndef SN_129
    #define SN_129
#endif // SN_129
#ifndef SB_129
    #define SB_129
#endif // SB_129
#ifndef TE_129
    #define TE_129
#endif // TE_129
#ifndef I_129
    #define I_129
#endif // I_129
#ifndef XE_129
    #define XE_129
#endif // XE_129
#ifndef CS_129
    #define CS_129
#endif // CS_129
#ifndef BA_129
    #define BA_129
#endif // BA_129
#ifndef LA_129
    #define LA_129
#endif // LA_129
#ifndef CE_129
    #define CE_129
#endif // CE_129
#ifndef PR_129
    #define PR_129
#endif // PR_129
#ifndef ND_129
    #define ND_129
#endif // ND_129
#ifndef PM_129
    #define PM_129
#endif // PM_129
#ifndef SM_129
    #define SM_129
#endif // SM_129
#ifndef PD_130
    #define PD_130
#endif // PD_130
#ifndef AG_130
    #define AG_130
#endif // AG_130
#ifndef CD_130
    #define CD_130
#endif // CD_130
#ifndef IN_130
    #define IN_130
#endif // IN_130
#ifndef SN_130
    #define SN_130
#endif // SN_130
#ifndef SB_130
    #define SB_130
#endif // SB_130
#ifndef TE_130
    #define TE_130
#endif // TE_130
#ifndef I_130
    #define I_130
#endif // I_130
#ifndef XE_130
    #define XE_130
#endif // XE_130
#ifndef CS_130
    #define CS_130
#endif // CS_130
#ifndef BA_130
    #define BA_130
#endif // BA_130
#ifndef LA_130
    #define LA_130
#endif // LA_130
#ifndef CE_130
    #define CE_130
#endif // CE_130
#ifndef PR_130
    #define PR_130
#endif // PR_130
#ifndef ND_130
    #define ND_130
#endif // ND_130
#ifndef PM_130
    #define PM_130
#endif // PM_130
#ifndef SM_130
    #define SM_130
#endif // SM_130
#ifndef EU_130
    #define EU_130
#endif // EU_130
#ifndef PD_131
    #define PD_131
#endif // PD_131
#ifndef AG_131
    #define AG_131
#endif // AG_131
#ifndef CD_131
    #define CD_131
#endif // CD_131
#ifndef IN_131
    #define IN_131
#endif // IN_131
#ifndef SN_131
    #define SN_131
#endif // SN_131
#ifndef SB_131
    #define SB_131
#endif // SB_131
#ifndef TE_131
    #define TE_131
#endif // TE_131
#ifndef I_131
    #define I_131
#endif // I_131
#ifndef XE_131
    #define XE_131
#endif // XE_131
#ifndef CS_131
    #define CS_131
#endif // CS_131
#ifndef BA_131
    #define BA_131
#endif // BA_131
#ifndef LA_131
    #define LA_131
#endif // LA_131
#ifndef CE_131
    #define CE_131
#endif // CE_131
#ifndef PR_131
    #define PR_131
#endif // PR_131
#ifndef ND_131
    #define ND_131
#endif // ND_131
#ifndef PM_131
    #define PM_131
#endif // PM_131
#ifndef SM_131
    #define SM_131
#endif // SM_131
#ifndef EU_131
    #define EU_131
#endif // EU_131
#ifndef AG_132
    #define AG_132
#endif // AG_132
#ifndef CD_132
    #define CD_132
#endif // CD_132
#ifndef IN_132
    #define IN_132
#endif // IN_132
#ifndef SN_132
    #define SN_132
#endif // SN_132
#ifndef SB_132
    #define SB_132
#endif // SB_132
#ifndef TE_132
    #define TE_132
#endif // TE_132
#ifndef I_132
    #define I_132
#endif // I_132
#ifndef XE_132
    #define XE_132
#endif // XE_132
#ifndef CS_132
    #define CS_132
#endif // CS_132
#ifndef BA_132
    #define BA_132
#endif // BA_132
#ifndef LA_132
    #define LA_132
#endif // LA_132
#ifndef CE_132
    #define CE_132
#endif // CE_132
#ifndef PR_132
    #define PR_132
#endif // PR_132
#ifndef ND_132
    #define ND_132
#endif // ND_132
#ifndef PM_132
    #define PM_132
#endif // PM_132
#ifndef SM_132
    #define SM_132
#endif // SM_132
#ifndef EU_132
    #define EU_132
#endif // EU_132
#ifndef AG_133
    #define AG_133
#endif // AG_133
#ifndef CD_133
    #define CD_133
#endif // CD_133
#ifndef IN_133
    #define IN_133
#endif // IN_133
#ifndef SN_133
    #define SN_133
#endif // SN_133
#ifndef SB_133
    #define SB_133
#endif // SB_133
#ifndef TE_133
    #define TE_133
#endif // TE_133
#ifndef I_133
    #define I_133
#endif // I_133
#ifndef XE_133
    #define XE_133
#endif // XE_133
#ifndef CS_133
    #define CS_133
#endif // CS_133
#ifndef BA_133
    #define BA_133
#endif // BA_133
#ifndef LA_133
    #define LA_133
#endif // LA_133
#ifndef CE_133
    #define CE_133
#endif // CE_133
#ifndef PR_133
    #define PR_133
#endif // PR_133
#ifndef ND_133
    #define ND_133
#endif // ND_133
#ifndef PM_133
    #define PM_133
#endif // PM_133
#ifndef SM_133
    #define SM_133
#endif // SM_133
#ifndef EU_133
    #define EU_133
#endif // EU_133
#ifndef GD_133
    #define GD_133
#endif // GD_133
#ifndef CD_134
    #define CD_134
#endif // CD_134
#ifndef IN_134
    #define IN_134
#endif // IN_134
#ifndef SN_134
    #define SN_134
#endif // SN_134
#ifndef SB_134
    #define SB_134
#endif // SB_134
#ifndef TE_134
    #define TE_134
#endif // TE_134
#ifndef I_134
    #define I_134
#endif // I_134
#ifndef XE_134
    #define XE_134
#endif // XE_134
#ifndef CS_134
    #define CS_134
#endif // CS_134
#ifndef BA_134
    #define BA_134
#endif // BA_134
#ifndef LA_134
    #define LA_134
#endif // LA_134
#ifndef CE_134
    #define CE_134
#endif // CE_134
#ifndef PR_134
    #define PR_134
#endif // PR_134
#ifndef ND_134
    #define ND_134
#endif // ND_134
#ifndef PM_134
    #define PM_134
#endif // PM_134
#ifndef SM_134
    #define SM_134
#endif // SM_134
#ifndef EU_134
    #define EU_134
#endif // EU_134
#ifndef GD_134
    #define GD_134
#endif // GD_134
#ifndef CD_135
    #define CD_135
#endif // CD_135
#ifndef IN_135
    #define IN_135
#endif // IN_135
#ifndef SN_135
    #define SN_135
#endif // SN_135
#ifndef SB_135
    #define SB_135
#endif // SB_135
#ifndef TE_135
    #define TE_135
#endif // TE_135
#ifndef I_135
    #define I_135
#endif // I_135
#ifndef XE_135
    #define XE_135
#endif // XE_135
#ifndef CS_135
    #define CS_135
#endif // CS_135
#ifndef BA_135
    #define BA_135
#endif // BA_135
#ifndef LA_135
    #define LA_135
#endif // LA_135
#ifndef CE_135
    #define CE_135
#endif // CE_135
#ifndef PR_135
    #define PR_135
#endif // PR_135
#ifndef ND_135
    #define ND_135
#endif // ND_135
#ifndef PM_135
    #define PM_135
#endif // PM_135
#ifndef SM_135
    #define SM_135
#endif // SM_135
#ifndef EU_135
    #define EU_135
#endif // EU_135
#ifndef GD_135
    #define GD_135
#endif // GD_135
#ifndef TB_135
    #define TB_135
#endif // TB_135
#ifndef IN_136
    #define IN_136
#endif // IN_136
#ifndef SN_136
    #define SN_136
#endif // SN_136
#ifndef SB_136
    #define SB_136
#endif // SB_136
#ifndef TE_136
    #define TE_136
#endif // TE_136
#ifndef I_136
    #define I_136
#endif // I_136
#ifndef XE_136
    #define XE_136
#endif // XE_136
#ifndef CS_136
    #define CS_136
#endif // CS_136
#ifndef BA_136
    #define BA_136
#endif // BA_136
#ifndef LA_136
    #define LA_136
#endif // LA_136
#ifndef CE_136
    #define CE_136
#endif // CE_136
#ifndef PR_136
    #define PR_136
#endif // PR_136
#ifndef ND_136
    #define ND_136
#endif // ND_136
#ifndef PM_136
    #define PM_136
#endif // PM_136
#ifndef SM_136
    #define SM_136
#endif // SM_136
#ifndef EU_136
    #define EU_136
#endif // EU_136
#ifndef GD_136
    #define GD_136
#endif // GD_136
#ifndef TB_136
    #define TB_136
#endif // TB_136
#ifndef IN_137
    #define IN_137
#endif // IN_137
#ifndef SN_137
    #define SN_137
#endif // SN_137
#ifndef SB_137
    #define SB_137
#endif // SB_137
#ifndef TE_137
    #define TE_137
#endif // TE_137
#ifndef I_137
    #define I_137
#endif // I_137
#ifndef XE_137
    #define XE_137
#endif // XE_137
#ifndef CS_137
    #define CS_137
#endif // CS_137
#ifndef BA_137
    #define BA_137
#endif // BA_137
#ifndef LA_137
    #define LA_137
#endif // LA_137
#ifndef CE_137
    #define CE_137
#endif // CE_137
#ifndef PR_137
    #define PR_137
#endif // PR_137
#ifndef ND_137
    #define ND_137
#endif // ND_137
#ifndef PM_137
    #define PM_137
#endif // PM_137
#ifndef SM_137
    #define SM_137
#endif // SM_137
#ifndef EU_137
    #define EU_137
#endif // EU_137
#ifndef GD_137
    #define GD_137
#endif // GD_137
#ifndef TB_137
    #define TB_137
#endif // TB_137
#ifndef SN_138
    #define SN_138
#endif // SN_138
#ifndef SB_138
    #define SB_138
#endif // SB_138
#ifndef TE_138
    #define TE_138
#endif // TE_138
#ifndef I_138
    #define I_138
#endif // I_138
#ifndef XE_138
    #define XE_138
#endif // XE_138
#ifndef CS_138
    #define CS_138
#endif // CS_138
#ifndef BA_138
    #define BA_138
#endif // BA_138
#ifndef LA_138
    #define LA_138
#endif // LA_138
#ifndef CE_138
    #define CE_138
#endif // CE_138
#ifndef PR_138
    #define PR_138
#endif // PR_138
#ifndef ND_138
    #define ND_138
#endif // ND_138
#ifndef PM_138
    #define PM_138
#endif // PM_138
#ifndef SM_138
    #define SM_138
#endif // SM_138
#ifndef EU_138
    #define EU_138
#endif // EU_138
#ifndef GD_138
    #define GD_138
#endif // GD_138
#ifndef TB_138
    #define TB_138
#endif // TB_138
#ifndef DY_138
    #define DY_138
#endif // DY_138
#ifndef SN_139
    #define SN_139
#endif // SN_139
#ifndef SB_139
    #define SB_139
#endif // SB_139
#ifndef TE_139
    #define TE_139
#endif // TE_139
#ifndef I_139
    #define I_139
#endif // I_139
#ifndef XE_139
    #define XE_139
#endif // XE_139
#ifndef CS_139
    #define CS_139
#endif // CS_139
#ifndef BA_139
    #define BA_139
#endif // BA_139
#ifndef LA_139
    #define LA_139
#endif // LA_139
#ifndef CE_139
    #define CE_139
#endif // CE_139
#ifndef PR_139
    #define PR_139
#endif // PR_139
#ifndef ND_139
    #define ND_139
#endif // ND_139
#ifndef PM_139
    #define PM_139
#endif // PM_139
#ifndef SM_139
    #define SM_139
#endif // SM_139
#ifndef EU_139
    #define EU_139
#endif // EU_139
#ifndef GD_139
    #define GD_139
#endif // GD_139
#ifndef TB_139
    #define TB_139
#endif // TB_139
#ifndef DY_139
    #define DY_139
#endif // DY_139
#ifndef SN_140
    #define SN_140
#endif // SN_140
#ifndef SB_140
    #define SB_140
#endif // SB_140
#ifndef TE_140
    #define TE_140
#endif // TE_140
#ifndef I_140
    #define I_140
#endif // I_140
#ifndef XE_140
    #define XE_140
#endif // XE_140
#ifndef CS_140
    #define CS_140
#endif // CS_140
#ifndef BA_140
    #define BA_140
#endif // BA_140
#ifndef LA_140
    #define LA_140
#endif // LA_140
#ifndef CE_140
    #define CE_140
#endif // CE_140
#ifndef PR_140
    #define PR_140
#endif // PR_140
#ifndef ND_140
    #define ND_140
#endif // ND_140
#ifndef PM_140
    #define PM_140
#endif // PM_140
#ifndef SM_140
    #define SM_140
#endif // SM_140
#ifndef EU_140
    #define EU_140
#endif // EU_140
#ifndef GD_140
    #define GD_140
#endif // GD_140
#ifndef TB_140
    #define TB_140
#endif // TB_140
#ifndef DY_140
    #define DY_140
#endif // DY_140
#ifndef HO_140
    #define HO_140
#endif // HO_140
#ifndef SB_141
    #define SB_141
#endif // SB_141
#ifndef TE_141
    #define TE_141
#endif // TE_141
#ifndef I_141
    #define I_141
#endif // I_141
#ifndef XE_141
    #define XE_141
#endif // XE_141
#ifndef CS_141
    #define CS_141
#endif // CS_141
#ifndef BA_141
    #define BA_141
#endif // BA_141
#ifndef LA_141
    #define LA_141
#endif // LA_141
#ifndef CE_141
    #define CE_141
#endif // CE_141
#ifndef PR_141
    #define PR_141
#endif // PR_141
#ifndef ND_141
    #define ND_141
#endif // ND_141
#ifndef PM_141
    #define PM_141
#endif // PM_141
#ifndef SM_141
    #define SM_141
#endif // SM_141
#ifndef EU_141
    #define EU_141
#endif // EU_141
#ifndef GD_141
    #define GD_141
#endif // GD_141
#ifndef TB_141
    #define TB_141
#endif // TB_141
#ifndef DY_141
    #define DY_141
#endif // DY_141
#ifndef HO_141
    #define HO_141
#endif // HO_141
#ifndef SB_142
    #define SB_142
#endif // SB_142
#ifndef TE_142
    #define TE_142
#endif // TE_142
#ifndef I_142
    #define I_142
#endif // I_142
#ifndef XE_142
    #define XE_142
#endif // XE_142
#ifndef CS_142
    #define CS_142
#endif // CS_142
#ifndef BA_142
    #define BA_142
#endif // BA_142
#ifndef LA_142
    #define LA_142
#endif // LA_142
#ifndef CE_142
    #define CE_142
#endif // CE_142
#ifndef PR_142
    #define PR_142
#endif // PR_142
#ifndef ND_142
    #define ND_142
#endif // ND_142
#ifndef PM_142
    #define PM_142
#endif // PM_142
#ifndef SM_142
    #define SM_142
#endif // SM_142
#ifndef EU_142
    #define EU_142
#endif // EU_142
#ifndef GD_142
    #define GD_142
#endif // GD_142
#ifndef TB_142
    #define TB_142
#endif // TB_142
#ifndef DY_142
    #define DY_142
#endif // DY_142
#ifndef HO_142
    #define HO_142
#endif // HO_142
#ifndef ER_142
    #define ER_142
#endif // ER_142
#ifndef TE_143
    #define TE_143
#endif // TE_143
#ifndef I_143
    #define I_143
#endif // I_143
#ifndef XE_143
    #define XE_143
#endif // XE_143
#ifndef CS_143
    #define CS_143
#endif // CS_143
#ifndef BA_143
    #define BA_143
#endif // BA_143
#ifndef LA_143
    #define LA_143
#endif // LA_143
#ifndef CE_143
    #define CE_143
#endif // CE_143
#ifndef PR_143
    #define PR_143
#endif // PR_143
#ifndef ND_143
    #define ND_143
#endif // ND_143
#ifndef PM_143
    #define PM_143
#endif // PM_143
#ifndef SM_143
    #define SM_143
#endif // SM_143
#ifndef EU_143
    #define EU_143
#endif // EU_143
#ifndef GD_143
    #define GD_143
#endif // GD_143
#ifndef TB_143
    #define TB_143
#endif // TB_143
#ifndef DY_143
    #define DY_143
#endif // DY_143
#ifndef HO_143
    #define HO_143
#endif // HO_143
#ifndef ER_143
    #define ER_143
#endif // ER_143
#ifndef TE_144
    #define TE_144
#endif // TE_144
#ifndef I_144
    #define I_144
#endif // I_144
#ifndef XE_144
    #define XE_144
#endif // XE_144
#ifndef CS_144
    #define CS_144
#endif // CS_144
#ifndef BA_144
    #define BA_144
#endif // BA_144
#ifndef LA_144
    #define LA_144
#endif // LA_144
#ifndef CE_144
    #define CE_144
#endif // CE_144
#ifndef PR_144
    #define PR_144
#endif // PR_144
#ifndef ND_144
    #define ND_144
#endif // ND_144
#ifndef PM_144
    #define PM_144
#endif // PM_144
#ifndef SM_144
    #define SM_144
#endif // SM_144
#ifndef EU_144
    #define EU_144
#endif // EU_144
#ifndef GD_144
    #define GD_144
#endif // GD_144
#ifndef TB_144
    #define TB_144
#endif // TB_144
#ifndef DY_144
    #define DY_144
#endif // DY_144
#ifndef HO_144
    #define HO_144
#endif // HO_144
#ifndef ER_144
    #define ER_144
#endif // ER_144
#ifndef TM_144
    #define TM_144
#endif // TM_144
#ifndef TE_145
    #define TE_145
#endif // TE_145
#ifndef I_145
    #define I_145
#endif // I_145
#ifndef XE_145
    #define XE_145
#endif // XE_145
#ifndef CS_145
    #define CS_145
#endif // CS_145
#ifndef BA_145
    #define BA_145
#endif // BA_145
#ifndef LA_145
    #define LA_145
#endif // LA_145
#ifndef CE_145
    #define CE_145
#endif // CE_145
#ifndef PR_145
    #define PR_145
#endif // PR_145
#ifndef ND_145
    #define ND_145
#endif // ND_145
#ifndef PM_145
    #define PM_145
#endif // PM_145
#ifndef SM_145
    #define SM_145
#endif // SM_145
#ifndef EU_145
    #define EU_145
#endif // EU_145
#ifndef GD_145
    #define GD_145
#endif // GD_145
#ifndef TB_145
    #define TB_145
#endif // TB_145
#ifndef DY_145
    #define DY_145
#endif // DY_145
#ifndef HO_145
    #define HO_145
#endif // HO_145
#ifndef ER_145
    #define ER_145
#endif // ER_145
#ifndef TM_145
    #define TM_145
#endif // TM_145
#ifndef I_146
    #define I_146
#endif // I_146
#ifndef XE_146
    #define XE_146
#endif // XE_146
#ifndef CS_146
    #define CS_146
#endif // CS_146
#ifndef BA_146
    #define BA_146
#endif // BA_146
#ifndef LA_146
    #define LA_146
#endif // LA_146
#ifndef CE_146
    #define CE_146
#endif // CE_146
#ifndef PR_146
    #define PR_146
#endif // PR_146
#ifndef ND_146
    #define ND_146
#endif // ND_146
#ifndef PM_146
    #define PM_146
#endif // PM_146
#ifndef SM_146
    #define SM_146
#endif // SM_146
#ifndef EU_146
    #define EU_146
#endif // EU_146
#ifndef GD_146
    #define GD_146
#endif // GD_146
#ifndef TB_146
    #define TB_146
#endif // TB_146
#ifndef DY_146
    #define DY_146
#endif // DY_146
#ifndef HO_146
    #define HO_146
#endif // HO_146
#ifndef ER_146
    #define ER_146
#endif // ER_146
#ifndef TM_146
    #define TM_146
#endif // TM_146
#ifndef I_147
    #define I_147
#endif // I_147
#ifndef XE_147
    #define XE_147
#endif // XE_147
#ifndef CS_147
    #define CS_147
#endif // CS_147
#ifndef BA_147
    #define BA_147
#endif // BA_147
#ifndef LA_147
    #define LA_147
#endif // LA_147
#ifndef CE_147
    #define CE_147
#endif // CE_147
#ifndef PR_147
    #define PR_147
#endif // PR_147
#ifndef ND_147
    #define ND_147
#endif // ND_147
#ifndef PM_147
    #define PM_147
#endif // PM_147
#ifndef SM_147
    #define SM_147
#endif // SM_147
#ifndef EU_147
    #define EU_147
#endif // EU_147
#ifndef GD_147
    #define GD_147
#endif // GD_147
#ifndef TB_147
    #define TB_147
#endif // TB_147
#ifndef DY_147
    #define DY_147
#endif // DY_147
#ifndef HO_147
    #define HO_147
#endif // HO_147
#ifndef ER_147
    #define ER_147
#endif // ER_147
#ifndef TM_147
    #define TM_147
#endif // TM_147
#ifndef XE_148
    #define XE_148
#endif // XE_148
#ifndef CS_148
    #define CS_148
#endif // CS_148
#ifndef BA_148
    #define BA_148
#endif // BA_148
#ifndef LA_148
    #define LA_148
#endif // LA_148
#ifndef CE_148
    #define CE_148
#endif // CE_148
#ifndef PR_148
    #define PR_148
#endif // PR_148
#ifndef ND_148
    #define ND_148
#endif // ND_148
#ifndef PM_148
    #define PM_148
#endif // PM_148
#ifndef SM_148
    #define SM_148
#endif // SM_148
#ifndef EU_148
    #define EU_148
#endif // EU_148
#ifndef GD_148
    #define GD_148
#endif // GD_148
#ifndef TB_148
    #define TB_148
#endif // TB_148
#ifndef DY_148
    #define DY_148
#endif // DY_148
#ifndef HO_148
    #define HO_148
#endif // HO_148
#ifndef ER_148
    #define ER_148
#endif // ER_148
#ifndef TM_148
    #define TM_148
#endif // TM_148
#ifndef YB_148
    #define YB_148
#endif // YB_148
#ifndef XE_149
    #define XE_149
#endif // XE_149
#ifndef CS_149
    #define CS_149
#endif // CS_149
#ifndef BA_149
    #define BA_149
#endif // BA_149
#ifndef LA_149
    #define LA_149
#endif // LA_149
#ifndef CE_149
    #define CE_149
#endif // CE_149
#ifndef PR_149
    #define PR_149
#endif // PR_149
#ifndef ND_149
    #define ND_149
#endif // ND_149
#ifndef PM_149
    #define PM_149
#endif // PM_149
#ifndef SM_149
    #define SM_149
#endif // SM_149
#ifndef EU_149
    #define EU_149
#endif // EU_149
#ifndef GD_149
    #define GD_149
#endif // GD_149
#ifndef TB_149
    #define TB_149
#endif // TB_149
#ifndef DY_149
    #define DY_149
#endif // DY_149
#ifndef HO_149
    #define HO_149
#endif // HO_149
#ifndef ER_149
    #define ER_149
#endif // ER_149
#ifndef TM_149
    #define TM_149
#endif // TM_149
#ifndef YB_149
    #define YB_149
#endif // YB_149
#ifndef XE_150
    #define XE_150
#endif // XE_150
#ifndef CS_150
    #define CS_150
#endif // CS_150
#ifndef BA_150
    #define BA_150
#endif // BA_150
#ifndef LA_150
    #define LA_150
#endif // LA_150
#ifndef CE_150
    #define CE_150
#endif // CE_150
#ifndef PR_150
    #define PR_150
#endif // PR_150
#ifndef ND_150
    #define ND_150
#endif // ND_150
#ifndef PM_150
    #define PM_150
#endif // PM_150
#ifndef SM_150
    #define SM_150
#endif // SM_150
#ifndef EU_150
    #define EU_150
#endif // EU_150
#ifndef GD_150
    #define GD_150
#endif // GD_150
#ifndef TB_150
    #define TB_150
#endif // TB_150
#ifndef DY_150
    #define DY_150
#endif // DY_150
#ifndef HO_150
    #define HO_150
#endif // HO_150
#ifndef ER_150
    #define ER_150
#endif // ER_150
#ifndef TM_150
    #define TM_150
#endif // TM_150
#ifndef YB_150
    #define YB_150
#endif // YB_150
#ifndef LU_150
    #define LU_150
#endif // LU_150
#ifndef CS_151
    #define CS_151
#endif // CS_151
#ifndef BA_151
    #define BA_151
#endif // BA_151
#ifndef LA_151
    #define LA_151
#endif // LA_151
#ifndef CE_151
    #define CE_151
#endif // CE_151
#ifndef PR_151
    #define PR_151
#endif // PR_151
#ifndef ND_151
    #define ND_151
#endif // ND_151
#ifndef PM_151
    #define PM_151
#endif // PM_151
#ifndef SM_151
    #define SM_151
#endif // SM_151
#ifndef EU_151
    #define EU_151
#endif // EU_151
#ifndef GD_151
    #define GD_151
#endif // GD_151
#ifndef TB_151
    #define TB_151
#endif // TB_151
#ifndef DY_151
    #define DY_151
#endif // DY_151
#ifndef HO_151
    #define HO_151
#endif // HO_151
#ifndef ER_151
    #define ER_151
#endif // ER_151
#ifndef TM_151
    #define TM_151
#endif // TM_151
#ifndef YB_151
    #define YB_151
#endif // YB_151
#ifndef LU_151
    #define LU_151
#endif // LU_151
#ifndef CS_152
    #define CS_152
#endif // CS_152
#ifndef BA_152
    #define BA_152
#endif // BA_152
#ifndef LA_152
    #define LA_152
#endif // LA_152
#ifndef CE_152
    #define CE_152
#endif // CE_152
#ifndef PR_152
    #define PR_152
#endif // PR_152
#ifndef ND_152
    #define ND_152
#endif // ND_152
#ifndef PM_152
    #define PM_152
#endif // PM_152
#ifndef SM_152
    #define SM_152
#endif // SM_152
#ifndef EU_152
    #define EU_152
#endif // EU_152
#ifndef GD_152
    #define GD_152
#endif // GD_152
#ifndef TB_152
    #define TB_152
#endif // TB_152
#ifndef DY_152
    #define DY_152
#endif // DY_152
#ifndef HO_152
    #define HO_152
#endif // HO_152
#ifndef ER_152
    #define ER_152
#endif // ER_152
#ifndef TM_152
    #define TM_152
#endif // TM_152
#ifndef YB_152
    #define YB_152
#endif // YB_152
#ifndef LU_152
    #define LU_152
#endif // LU_152
#ifndef BA_153
    #define BA_153
#endif // BA_153
#ifndef LA_153
    #define LA_153
#endif // LA_153
#ifndef CE_153
    #define CE_153
#endif // CE_153
#ifndef PR_153
    #define PR_153
#endif // PR_153
#ifndef ND_153
    #define ND_153
#endif // ND_153
#ifndef PM_153
    #define PM_153
#endif // PM_153
#ifndef SM_153
    #define SM_153
#endif // SM_153
#ifndef EU_153
    #define EU_153
#endif // EU_153
#ifndef GD_153
    #define GD_153
#endif // GD_153
#ifndef TB_153
    #define TB_153
#endif // TB_153
#ifndef DY_153
    #define DY_153
#endif // DY_153
#ifndef HO_153
    #define HO_153
#endif // HO_153
#ifndef ER_153
    #define ER_153
#endif // ER_153
#ifndef TM_153
    #define TM_153
#endif // TM_153
#ifndef YB_153
    #define YB_153
#endif // YB_153
#ifndef LU_153
    #define LU_153
#endif // LU_153
#ifndef HF_153
    #define HF_153
#endif // HF_153
#ifndef BA_154
    #define BA_154
#endif // BA_154
#ifndef LA_154
    #define LA_154
#endif // LA_154
#ifndef CE_154
    #define CE_154
#endif // CE_154
#ifndef PR_154
    #define PR_154
#endif // PR_154
#ifndef ND_154
    #define ND_154
#endif // ND_154
#ifndef PM_154
    #define PM_154
#endif // PM_154
#ifndef SM_154
    #define SM_154
#endif // SM_154
#ifndef EU_154
    #define EU_154
#endif // EU_154
#ifndef GD_154
    #define GD_154
#endif // GD_154
#ifndef TB_154
    #define TB_154
#endif // TB_154
#ifndef DY_154
    #define DY_154
#endif // DY_154
#ifndef HO_154
    #define HO_154
#endif // HO_154
#ifndef ER_154
    #define ER_154
#endif // ER_154
#ifndef TM_154
    #define TM_154
#endif // TM_154
#ifndef YB_154
    #define YB_154
#endif // YB_154
#ifndef LU_154
    #define LU_154
#endif // LU_154
#ifndef HF_154
    #define HF_154
#endif // HF_154
#ifndef LA_155
    #define LA_155
#endif // LA_155
#ifndef CE_155
    #define CE_155
#endif // CE_155
#ifndef PR_155
    #define PR_155
#endif // PR_155
#ifndef ND_155
    #define ND_155
#endif // ND_155
#ifndef PM_155
    #define PM_155
#endif // PM_155
#ifndef SM_155
    #define SM_155
#endif // SM_155
#ifndef EU_155
    #define EU_155
#endif // EU_155
#ifndef GD_155
    #define GD_155
#endif // GD_155
#ifndef TB_155
    #define TB_155
#endif // TB_155
#ifndef DY_155
    #define DY_155
#endif // DY_155
#ifndef HO_155
    #define HO_155
#endif // HO_155
#ifndef ER_155
    #define ER_155
#endif // ER_155
#ifndef TM_155
    #define TM_155
#endif // TM_155
#ifndef YB_155
    #define YB_155
#endif // YB_155
#ifndef LU_155
    #define LU_155
#endif // LU_155
#ifndef HF_155
    #define HF_155
#endif // HF_155
#ifndef TA_155
    #define TA_155
#endif // TA_155
#ifndef LA_156
    #define LA_156
#endif // LA_156
#ifndef CE_156
    #define CE_156
#endif // CE_156
#ifndef PR_156
    #define PR_156
#endif // PR_156
#ifndef ND_156
    #define ND_156
#endif // ND_156
#ifndef PM_156
    #define PM_156
#endif // PM_156
#ifndef SM_156
    #define SM_156
#endif // SM_156
#ifndef EU_156
    #define EU_156
#endif // EU_156
#ifndef GD_156
    #define GD_156
#endif // GD_156
#ifndef TB_156
    #define TB_156
#endif // TB_156
#ifndef DY_156
    #define DY_156
#endif // DY_156
#ifndef HO_156
    #define HO_156
#endif // HO_156
#ifndef ER_156
    #define ER_156
#endif // ER_156
#ifndef TM_156
    #define TM_156
#endif // TM_156
#ifndef YB_156
    #define YB_156
#endif // YB_156
#ifndef LU_156
    #define LU_156
#endif // LU_156
#ifndef HF_156
    #define HF_156
#endif // HF_156
#ifndef TA_156
    #define TA_156
#endif // TA_156
#ifndef LA_157
    #define LA_157
#endif // LA_157
#ifndef CE_157
    #define CE_157
#endif // CE_157
#ifndef PR_157
    #define PR_157
#endif // PR_157
#ifndef ND_157
    #define ND_157
#endif // ND_157
#ifndef PM_157
    #define PM_157
#endif // PM_157
#ifndef SM_157
    #define SM_157
#endif // SM_157
#ifndef EU_157
    #define EU_157
#endif // EU_157
#ifndef GD_157
    #define GD_157
#endif // GD_157
#ifndef TB_157
    #define TB_157
#endif // TB_157
#ifndef DY_157
    #define DY_157
#endif // DY_157
#ifndef HO_157
    #define HO_157
#endif // HO_157
#ifndef ER_157
    #define ER_157
#endif // ER_157
#ifndef TM_157
    #define TM_157
#endif // TM_157
#ifndef YB_157
    #define YB_157
#endif // YB_157
#ifndef LU_157
    #define LU_157
#endif // LU_157
#ifndef HF_157
    #define HF_157
#endif // HF_157
#ifndef TA_157
    #define TA_157
#endif // TA_157
#ifndef W_157
    #define W_157
#endif // W_157
#ifndef CE_158
    #define CE_158
#endif // CE_158
#ifndef PR_158
    #define PR_158
#endif // PR_158
#ifndef ND_158
    #define ND_158
#endif // ND_158
#ifndef PM_158
    #define PM_158
#endif // PM_158
#ifndef SM_158
    #define SM_158
#endif // SM_158
#ifndef EU_158
    #define EU_158
#endif // EU_158
#ifndef GD_158
    #define GD_158
#endif // GD_158
#ifndef TB_158
    #define TB_158
#endif // TB_158
#ifndef DY_158
    #define DY_158
#endif // DY_158
#ifndef HO_158
    #define HO_158
#endif // HO_158
#ifndef ER_158
    #define ER_158
#endif // ER_158
#ifndef TM_158
    #define TM_158
#endif // TM_158
#ifndef YB_158
    #define YB_158
#endif // YB_158
#ifndef LU_158
    #define LU_158
#endif // LU_158
#ifndef HF_158
    #define HF_158
#endif // HF_158
#ifndef TA_158
    #define TA_158
#endif // TA_158
#ifndef W_158
    #define W_158
#endif // W_158
#ifndef CE_159
    #define CE_159
#endif // CE_159
#ifndef PR_159
    #define PR_159
#endif // PR_159
#ifndef ND_159
    #define ND_159
#endif // ND_159
#ifndef PM_159
    #define PM_159
#endif // PM_159
#ifndef SM_159
    #define SM_159
#endif // SM_159
#ifndef EU_159
    #define EU_159
#endif // EU_159
#ifndef GD_159
    #define GD_159
#endif // GD_159
#ifndef TB_159
    #define TB_159
#endif // TB_159
#ifndef DY_159
    #define DY_159
#endif // DY_159
#ifndef HO_159
    #define HO_159
#endif // HO_159
#ifndef ER_159
    #define ER_159
#endif // ER_159
#ifndef TM_159
    #define TM_159
#endif // TM_159
#ifndef YB_159
    #define YB_159
#endif // YB_159
#ifndef LU_159
    #define LU_159
#endif // LU_159
#ifndef HF_159
    #define HF_159
#endif // HF_159
#ifndef TA_159
    #define TA_159
#endif // TA_159
#ifndef W_159
    #define W_159
#endif // W_159
#ifndef RE_159
    #define RE_159
#endif // RE_159
#ifndef PR_160
    #define PR_160
#endif // PR_160
#ifndef ND_160
    #define ND_160
#endif // ND_160
#ifndef PM_160
    #define PM_160
#endif // PM_160
#ifndef SM_160
    #define SM_160
#endif // SM_160
#ifndef EU_160
    #define EU_160
#endif // EU_160
#ifndef GD_160
    #define GD_160
#endif // GD_160
#ifndef TB_160
    #define TB_160
#endif // TB_160
#ifndef DY_160
    #define DY_160
#endif // DY_160
#ifndef HO_160
    #define HO_160
#endif // HO_160
#ifndef ER_160
    #define ER_160
#endif // ER_160
#ifndef TM_160
    #define TM_160
#endif // TM_160
#ifndef YB_160
    #define YB_160
#endif // YB_160
#ifndef LU_160
    #define LU_160
#endif // LU_160
#ifndef HF_160
    #define HF_160
#endif // HF_160
#ifndef TA_160
    #define TA_160
#endif // TA_160
#ifndef W_160
    #define W_160
#endif // W_160
#ifndef RE_160
    #define RE_160
#endif // RE_160
#ifndef PR_161
    #define PR_161
#endif // PR_161
#ifndef ND_161
    #define ND_161
#endif // ND_161
#ifndef PM_161
    #define PM_161
#endif // PM_161
#ifndef SM_161
    #define SM_161
#endif // SM_161
#ifndef EU_161
    #define EU_161
#endif // EU_161
#ifndef GD_161
    #define GD_161
#endif // GD_161
#ifndef TB_161
    #define TB_161
#endif // TB_161
#ifndef DY_161
    #define DY_161
#endif // DY_161
#ifndef HO_161
    #define HO_161
#endif // HO_161
#ifndef ER_161
    #define ER_161
#endif // ER_161
#ifndef TM_161
    #define TM_161
#endif // TM_161
#ifndef YB_161
    #define YB_161
#endif // YB_161
#ifndef LU_161
    #define LU_161
#endif // LU_161
#ifndef HF_161
    #define HF_161
#endif // HF_161
#ifndef TA_161
    #define TA_161
#endif // TA_161
#ifndef W_161
    #define W_161
#endif // W_161
#ifndef RE_161
    #define RE_161
#endif // RE_161
#ifndef OS_161
    #define OS_161
#endif // OS_161
#ifndef ND_162
    #define ND_162
#endif // ND_162
#ifndef PM_162
    #define PM_162
#endif // PM_162
#ifndef SM_162
    #define SM_162
#endif // SM_162
#ifndef EU_162
    #define EU_162
#endif // EU_162
#ifndef GD_162
    #define GD_162
#endif // GD_162
#ifndef TB_162
    #define TB_162
#endif // TB_162
#ifndef DY_162
    #define DY_162
#endif // DY_162
#ifndef HO_162
    #define HO_162
#endif // HO_162
#ifndef ER_162
    #define ER_162
#endif // ER_162
#ifndef TM_162
    #define TM_162
#endif // TM_162
#ifndef YB_162
    #define YB_162
#endif // YB_162
#ifndef LU_162
    #define LU_162
#endif // LU_162
#ifndef HF_162
    #define HF_162
#endif // HF_162
#ifndef TA_162
    #define TA_162
#endif // TA_162
#ifndef W_162
    #define W_162
#endif // W_162
#ifndef RE_162
    #define RE_162
#endif // RE_162
#ifndef OS_162
    #define OS_162
#endif // OS_162
#ifndef ND_163
    #define ND_163
#endif // ND_163
#ifndef PM_163
    #define PM_163
#endif // PM_163
#ifndef SM_163
    #define SM_163
#endif // SM_163
#ifndef EU_163
    #define EU_163
#endif // EU_163
#ifndef GD_163
    #define GD_163
#endif // GD_163
#ifndef TB_163
    #define TB_163
#endif // TB_163
#ifndef DY_163
    #define DY_163
#endif // DY_163
#ifndef HO_163
    #define HO_163
#endif // HO_163
#ifndef ER_163
    #define ER_163
#endif // ER_163
#ifndef TM_163
    #define TM_163
#endif // TM_163
#ifndef YB_163
    #define YB_163
#endif // YB_163
#ifndef LU_163
    #define LU_163
#endif // LU_163
#ifndef HF_163
    #define HF_163
#endif // HF_163
#ifndef TA_163
    #define TA_163
#endif // TA_163
#ifndef W_163
    #define W_163
#endif // W_163
#ifndef RE_163
    #define RE_163
#endif // RE_163
#ifndef OS_163
    #define OS_163
#endif // OS_163
#ifndef IR_163
    #define IR_163
#endif // IR_163
#ifndef PM_164
    #define PM_164
#endif // PM_164
#ifndef SM_164
    #define SM_164
#endif // SM_164
#ifndef EU_164
    #define EU_164
#endif // EU_164
#ifndef GD_164
    #define GD_164
#endif // GD_164
#ifndef TB_164
    #define TB_164
#endif // TB_164
#ifndef DY_164
    #define DY_164
#endif // DY_164
#ifndef HO_164
    #define HO_164
#endif // HO_164
#ifndef ER_164
    #define ER_164
#endif // ER_164
#ifndef TM_164
    #define TM_164
#endif // TM_164
#ifndef YB_164
    #define YB_164
#endif // YB_164
#ifndef LU_164
    #define LU_164
#endif // LU_164
#ifndef HF_164
    #define HF_164
#endif // HF_164
#ifndef TA_164
    #define TA_164
#endif // TA_164
#ifndef W_164
    #define W_164
#endif // W_164
#ifndef RE_164
    #define RE_164
#endif // RE_164
#ifndef OS_164
    #define OS_164
#endif // OS_164
#ifndef IR_164
    #define IR_164
#endif // IR_164
#ifndef PM_165
    #define PM_165
#endif // PM_165
#ifndef SM_165
    #define SM_165
#endif // SM_165
#ifndef EU_165
    #define EU_165
#endif // EU_165
#ifndef GD_165
    #define GD_165
#endif // GD_165
#ifndef TB_165
    #define TB_165
#endif // TB_165
#ifndef DY_165
    #define DY_165
#endif // DY_165
#ifndef HO_165
    #define HO_165
#endif // HO_165
#ifndef ER_165
    #define ER_165
#endif // ER_165
#ifndef TM_165
    #define TM_165
#endif // TM_165
#ifndef YB_165
    #define YB_165
#endif // YB_165
#ifndef LU_165
    #define LU_165
#endif // LU_165
#ifndef HF_165
    #define HF_165
#endif // HF_165
#ifndef TA_165
    #define TA_165
#endif // TA_165
#ifndef W_165
    #define W_165
#endif // W_165
#ifndef RE_165
    #define RE_165
#endif // RE_165
#ifndef OS_165
    #define OS_165
#endif // OS_165
#ifndef IR_165
    #define IR_165
#endif // IR_165
#ifndef PT_165
    #define PT_165
#endif // PT_165
#ifndef SM_166
    #define SM_166
#endif // SM_166
#ifndef EU_166
    #define EU_166
#endif // EU_166
#ifndef GD_166
    #define GD_166
#endif // GD_166
#ifndef TB_166
    #define TB_166
#endif // TB_166
#ifndef DY_166
    #define DY_166
#endif // DY_166
#ifndef HO_166
    #define HO_166
#endif // HO_166
#ifndef ER_166
    #define ER_166
#endif // ER_166
#ifndef TM_166
    #define TM_166
#endif // TM_166
#ifndef YB_166
    #define YB_166
#endif // YB_166
#ifndef LU_166
    #define LU_166
#endif // LU_166
#ifndef HF_166
    #define HF_166
#endif // HF_166
#ifndef TA_166
    #define TA_166
#endif // TA_166
#ifndef W_166
    #define W_166
#endif // W_166
#ifndef RE_166
    #define RE_166
#endif // RE_166
#ifndef OS_166
    #define OS_166
#endif // OS_166
#ifndef IR_166
    #define IR_166
#endif // IR_166
#ifndef PT_166
    #define PT_166
#endif // PT_166
#ifndef SM_167
    #define SM_167
#endif // SM_167
#ifndef EU_167
    #define EU_167
#endif // EU_167
#ifndef GD_167
    #define GD_167
#endif // GD_167
#ifndef TB_167
    #define TB_167
#endif // TB_167
#ifndef DY_167
    #define DY_167
#endif // DY_167
#ifndef HO_167
    #define HO_167
#endif // HO_167
#ifndef ER_167
    #define ER_167
#endif // ER_167
#ifndef TM_167
    #define TM_167
#endif // TM_167
#ifndef YB_167
    #define YB_167
#endif // YB_167
#ifndef LU_167
    #define LU_167
#endif // LU_167
#ifndef HF_167
    #define HF_167
#endif // HF_167
#ifndef TA_167
    #define TA_167
#endif // TA_167
#ifndef W_167
    #define W_167
#endif // W_167
#ifndef RE_167
    #define RE_167
#endif // RE_167
#ifndef OS_167
    #define OS_167
#endif // OS_167
#ifndef IR_167
    #define IR_167
#endif // IR_167
#ifndef PT_167
    #define PT_167
#endif // PT_167
#ifndef SM_168
    #define SM_168
#endif // SM_168
#ifndef EU_168
    #define EU_168
#endif // EU_168
#ifndef GD_168
    #define GD_168
#endif // GD_168
#ifndef TB_168
    #define TB_168
#endif // TB_168
#ifndef DY_168
    #define DY_168
#endif // DY_168
#ifndef HO_168
    #define HO_168
#endif // HO_168
#ifndef ER_168
    #define ER_168
#endif // ER_168
#ifndef TM_168
    #define TM_168
#endif // TM_168
#ifndef YB_168
    #define YB_168
#endif // YB_168
#ifndef LU_168
    #define LU_168
#endif // LU_168
#ifndef HF_168
    #define HF_168
#endif // HF_168
#ifndef TA_168
    #define TA_168
#endif // TA_168
#ifndef W_168
    #define W_168
#endif // W_168
#ifndef RE_168
    #define RE_168
#endif // RE_168
#ifndef OS_168
    #define OS_168
#endif // OS_168
#ifndef IR_168
    #define IR_168
#endif // IR_168
#ifndef PT_168
    #define PT_168
#endif // PT_168
#ifndef AU_168
    #define AU_168
#endif // AU_168
#ifndef EU_169
    #define EU_169
#endif // EU_169
#ifndef GD_169
    #define GD_169
#endif // GD_169
#ifndef TB_169
    #define TB_169
#endif // TB_169
#ifndef DY_169
    #define DY_169
#endif // DY_169
#ifndef HO_169
    #define HO_169
#endif // HO_169
#ifndef ER_169
    #define ER_169
#endif // ER_169
#ifndef TM_169
    #define TM_169
#endif // TM_169
#ifndef YB_169
    #define YB_169
#endif // YB_169
#ifndef LU_169
    #define LU_169
#endif // LU_169
#ifndef HF_169
    #define HF_169
#endif // HF_169
#ifndef TA_169
    #define TA_169
#endif // TA_169
#ifndef W_169
    #define W_169
#endif // W_169
#ifndef RE_169
    #define RE_169
#endif // RE_169
#ifndef OS_169
    #define OS_169
#endif // OS_169
#ifndef IR_169
    #define IR_169
#endif // IR_169
#ifndef PT_169
    #define PT_169
#endif // PT_169
#ifndef AU_169
    #define AU_169
#endif // AU_169
#ifndef EU_170
    #define EU_170
#endif // EU_170
#ifndef GD_170
    #define GD_170
#endif // GD_170
#ifndef TB_170
    #define TB_170
#endif // TB_170
#ifndef DY_170
    #define DY_170
#endif // DY_170
#ifndef HO_170
    #define HO_170
#endif // HO_170
#ifndef ER_170
    #define ER_170
#endif // ER_170
#ifndef TM_170
    #define TM_170
#endif // TM_170
#ifndef YB_170
    #define YB_170
#endif // YB_170
#ifndef LU_170
    #define LU_170
#endif // LU_170
#ifndef HF_170
    #define HF_170
#endif // HF_170
#ifndef TA_170
    #define TA_170
#endif // TA_170
#ifndef W_170
    #define W_170
#endif // W_170
#ifndef RE_170
    #define RE_170
#endif // RE_170
#ifndef OS_170
    #define OS_170
#endif // OS_170
#ifndef IR_170
    #define IR_170
#endif // IR_170
#ifndef PT_170
    #define PT_170
#endif // PT_170
#ifndef AU_170
    #define AU_170
#endif // AU_170
#ifndef HG_170
    #define HG_170
#endif // HG_170
#ifndef GD_171
    #define GD_171
#endif // GD_171
#ifndef TB_171
    #define TB_171
#endif // TB_171
#ifndef DY_171
    #define DY_171
#endif // DY_171
#ifndef HO_171
    #define HO_171
#endif // HO_171
#ifndef ER_171
    #define ER_171
#endif // ER_171
#ifndef TM_171
    #define TM_171
#endif // TM_171
#ifndef YB_171
    #define YB_171
#endif // YB_171
#ifndef LU_171
    #define LU_171
#endif // LU_171
#ifndef HF_171
    #define HF_171
#endif // HF_171
#ifndef TA_171
    #define TA_171
#endif // TA_171
#ifndef W_171
    #define W_171
#endif // W_171
#ifndef RE_171
    #define RE_171
#endif // RE_171
#ifndef OS_171
    #define OS_171
#endif // OS_171
#ifndef IR_171
    #define IR_171
#endif // IR_171
#ifndef PT_171
    #define PT_171
#endif // PT_171
#ifndef AU_171
    #define AU_171
#endif // AU_171
#ifndef HG_171
    #define HG_171
#endif // HG_171
#ifndef GD_172
    #define GD_172
#endif // GD_172
#ifndef TB_172
    #define TB_172
#endif // TB_172
#ifndef DY_172
    #define DY_172
#endif // DY_172
#ifndef HO_172
    #define HO_172
#endif // HO_172
#ifndef ER_172
    #define ER_172
#endif // ER_172
#ifndef TM_172
    #define TM_172
#endif // TM_172
#ifndef YB_172
    #define YB_172
#endif // YB_172
#ifndef LU_172
    #define LU_172
#endif // LU_172
#ifndef HF_172
    #define HF_172
#endif // HF_172
#ifndef TA_172
    #define TA_172
#endif // TA_172
#ifndef W_172
    #define W_172
#endif // W_172
#ifndef RE_172
    #define RE_172
#endif // RE_172
#ifndef OS_172
    #define OS_172
#endif // OS_172
#ifndef IR_172
    #define IR_172
#endif // IR_172
#ifndef PT_172
    #define PT_172
#endif // PT_172
#ifndef AU_172
    #define AU_172
#endif // AU_172
#ifndef HG_172
    #define HG_172
#endif // HG_172
#ifndef TB_173
    #define TB_173
#endif // TB_173
#ifndef DY_173
    #define DY_173
#endif // DY_173
#ifndef HO_173
    #define HO_173
#endif // HO_173
#ifndef ER_173
    #define ER_173
#endif // ER_173
#ifndef TM_173
    #define TM_173
#endif // TM_173
#ifndef YB_173
    #define YB_173
#endif // YB_173
#ifndef LU_173
    #define LU_173
#endif // LU_173
#ifndef HF_173
    #define HF_173
#endif // HF_173
#ifndef TA_173
    #define TA_173
#endif // TA_173
#ifndef W_173
    #define W_173
#endif // W_173
#ifndef RE_173
    #define RE_173
#endif // RE_173
#ifndef OS_173
    #define OS_173
#endif // OS_173
#ifndef IR_173
    #define IR_173
#endif // IR_173
#ifndef PT_173
    #define PT_173
#endif // PT_173
#ifndef AU_173
    #define AU_173
#endif // AU_173
#ifndef HG_173
    #define HG_173
#endif // HG_173
#ifndef TB_174
    #define TB_174
#endif // TB_174
#ifndef DY_174
    #define DY_174
#endif // DY_174
#ifndef HO_174
    #define HO_174
#endif // HO_174
#ifndef ER_174
    #define ER_174
#endif // ER_174
#ifndef TM_174
    #define TM_174
#endif // TM_174
#ifndef YB_174
    #define YB_174
#endif // YB_174
#ifndef LU_174
    #define LU_174
#endif // LU_174
#ifndef HF_174
    #define HF_174
#endif // HF_174
#ifndef TA_174
    #define TA_174
#endif // TA_174
#ifndef W_174
    #define W_174
#endif // W_174
#ifndef RE_174
    #define RE_174
#endif // RE_174
#ifndef OS_174
    #define OS_174
#endif // OS_174
#ifndef IR_174
    #define IR_174
#endif // IR_174
#ifndef PT_174
    #define PT_174
#endif // PT_174
#ifndef AU_174
    #define AU_174
#endif // AU_174
#ifndef HG_174
    #define HG_174
#endif // HG_174
#ifndef DY_175
    #define DY_175
#endif // DY_175
#ifndef HO_175
    #define HO_175
#endif // HO_175
#ifndef ER_175
    #define ER_175
#endif // ER_175
#ifndef TM_175
    #define TM_175
#endif // TM_175
#ifndef YB_175
    #define YB_175
#endif // YB_175
#ifndef LU_175
    #define LU_175
#endif // LU_175
#ifndef HF_175
    #define HF_175
#endif // HF_175
#ifndef TA_175
    #define TA_175
#endif // TA_175
#ifndef W_175
    #define W_175
#endif // W_175
#ifndef RE_175
    #define RE_175
#endif // RE_175
#ifndef OS_175
    #define OS_175
#endif // OS_175
#ifndef IR_175
    #define IR_175
#endif // IR_175
#ifndef PT_175
    #define PT_175
#endif // PT_175
#ifndef AU_175
    #define AU_175
#endif // AU_175
#ifndef HG_175
    #define HG_175
#endif // HG_175
#ifndef DY_176
    #define DY_176
#endif // DY_176
#ifndef HO_176
    #define HO_176
#endif // HO_176
#ifndef ER_176
    #define ER_176
#endif // ER_176
#ifndef TM_176
    #define TM_176
#endif // TM_176
#ifndef YB_176
    #define YB_176
#endif // YB_176
#ifndef LU_176
    #define LU_176
#endif // LU_176
#ifndef HF_176
    #define HF_176
#endif // HF_176
#ifndef TA_176
    #define TA_176
#endif // TA_176
#ifndef W_176
    #define W_176
#endif // W_176
#ifndef RE_176
    #define RE_176
#endif // RE_176
#ifndef OS_176
    #define OS_176
#endif // OS_176
#ifndef IR_176
    #define IR_176
#endif // IR_176
#ifndef PT_176
    #define PT_176
#endif // PT_176
#ifndef AU_176
    #define AU_176
#endif // AU_176
#ifndef HG_176
    #define HG_176
#endif // HG_176
#ifndef TL_176
    #define TL_176
#endif // TL_176
#ifndef HO_177
    #define HO_177
#endif // HO_177
#ifndef ER_177
    #define ER_177
#endif // ER_177
#ifndef TM_177
    #define TM_177
#endif // TM_177
#ifndef YB_177
    #define YB_177
#endif // YB_177
#ifndef LU_177
    #define LU_177
#endif // LU_177
#ifndef HF_177
    #define HF_177
#endif // HF_177
#ifndef TA_177
    #define TA_177
#endif // TA_177
#ifndef W_177
    #define W_177
#endif // W_177
#ifndef RE_177
    #define RE_177
#endif // RE_177
#ifndef OS_177
    #define OS_177
#endif // OS_177
#ifndef IR_177
    #define IR_177
#endif // IR_177
#ifndef PT_177
    #define PT_177
#endif // PT_177
#ifndef AU_177
    #define AU_177
#endif // AU_177
#ifndef HG_177
    #define HG_177
#endif // HG_177
#ifndef TL_177
    #define TL_177
#endif // TL_177
#ifndef HO_178
    #define HO_178
#endif // HO_178
#ifndef ER_178
    #define ER_178
#endif // ER_178
#ifndef TM_178
    #define TM_178
#endif // TM_178
#ifndef YB_178
    #define YB_178
#endif // YB_178
#ifndef LU_178
    #define LU_178
#endif // LU_178
#ifndef HF_178
    #define HF_178
#endif // HF_178
#ifndef TA_178
    #define TA_178
#endif // TA_178
#ifndef W_178
    #define W_178
#endif // W_178
#ifndef RE_178
    #define RE_178
#endif // RE_178
#ifndef OS_178
    #define OS_178
#endif // OS_178
#ifndef IR_178
    #define IR_178
#endif // IR_178
#ifndef PT_178
    #define PT_178
#endif // PT_178
#ifndef AU_178
    #define AU_178
#endif // AU_178
#ifndef HG_178
    #define HG_178
#endif // HG_178
#ifndef TL_178
    #define TL_178
#endif // TL_178
#ifndef PB_178
    #define PB_178
#endif // PB_178
#ifndef ER_179
    #define ER_179
#endif // ER_179
#ifndef TM_179
    #define TM_179
#endif // TM_179
#ifndef YB_179
    #define YB_179
#endif // YB_179
#ifndef LU_179
    #define LU_179
#endif // LU_179
#ifndef HF_179
    #define HF_179
#endif // HF_179
#ifndef TA_179
    #define TA_179
#endif // TA_179
#ifndef W_179
    #define W_179
#endif // W_179
#ifndef RE_179
    #define RE_179
#endif // RE_179
#ifndef OS_179
    #define OS_179
#endif // OS_179
#ifndef IR_179
    #define IR_179
#endif // IR_179
#ifndef PT_179
    #define PT_179
#endif // PT_179
#ifndef AU_179
    #define AU_179
#endif // AU_179
#ifndef HG_179
    #define HG_179
#endif // HG_179
#ifndef TL_179
    #define TL_179
#endif // TL_179
#ifndef PB_179
    #define PB_179
#endif // PB_179
#ifndef ER_180
    #define ER_180
#endif // ER_180
#ifndef TM_180
    #define TM_180
#endif // TM_180
#ifndef YB_180
    #define YB_180
#endif // YB_180
#ifndef LU_180
    #define LU_180
#endif // LU_180
#ifndef HF_180
    #define HF_180
#endif // HF_180
#ifndef TA_180
    #define TA_180
#endif // TA_180
#ifndef W_180
    #define W_180
#endif // W_180
#ifndef RE_180
    #define RE_180
#endif // RE_180
#ifndef OS_180
    #define OS_180
#endif // OS_180
#ifndef IR_180
    #define IR_180
#endif // IR_180
#ifndef PT_180
    #define PT_180
#endif // PT_180
#ifndef AU_180
    #define AU_180
#endif // AU_180
#ifndef HG_180
    #define HG_180
#endif // HG_180
#ifndef TL_180
    #define TL_180
#endif // TL_180
#ifndef PB_180
    #define PB_180
#endif // PB_180
#ifndef TM_181
    #define TM_181
#endif // TM_181
#ifndef YB_181
    #define YB_181
#endif // YB_181
#ifndef LU_181
    #define LU_181
#endif // LU_181
#ifndef HF_181
    #define HF_181
#endif // HF_181
#ifndef TA_181
    #define TA_181
#endif // TA_181
#ifndef W_181
    #define W_181
#endif // W_181
#ifndef RE_181
    #define RE_181
#endif // RE_181
#ifndef OS_181
    #define OS_181
#endif // OS_181
#ifndef IR_181
    #define IR_181
#endif // IR_181
#ifndef PT_181
    #define PT_181
#endif // PT_181
#ifndef AU_181
    #define AU_181
#endif // AU_181
#ifndef HG_181
    #define HG_181
#endif // HG_181
#ifndef TL_181
    #define TL_181
#endif // TL_181
#ifndef PB_181
    #define PB_181
#endif // PB_181
#ifndef TM_182
    #define TM_182
#endif // TM_182
#ifndef YB_182
    #define YB_182
#endif // YB_182
#ifndef LU_182
    #define LU_182
#endif // LU_182
#ifndef HF_182
    #define HF_182
#endif // HF_182
#ifndef TA_182
    #define TA_182
#endif // TA_182
#ifndef W_182
    #define W_182
#endif // W_182
#ifndef RE_182
    #define RE_182
#endif // RE_182
#ifndef OS_182
    #define OS_182
#endif // OS_182
#ifndef IR_182
    #define IR_182
#endif // IR_182
#ifndef PT_182
    #define PT_182
#endif // PT_182
#ifndef AU_182
    #define AU_182
#endif // AU_182
#ifndef HG_182
    #define HG_182
#endif // HG_182
#ifndef TL_182
    #define TL_182
#endif // TL_182
#ifndef PB_182
    #define PB_182
#endif // PB_182
#ifndef YB_183
    #define YB_183
#endif // YB_183
#ifndef LU_183
    #define LU_183
#endif // LU_183
#ifndef HF_183
    #define HF_183
#endif // HF_183
#ifndef TA_183
    #define TA_183
#endif // TA_183
#ifndef W_183
    #define W_183
#endif // W_183
#ifndef RE_183
    #define RE_183
#endif // RE_183
#ifndef OS_183
    #define OS_183
#endif // OS_183
#ifndef IR_183
    #define IR_183
#endif // IR_183
#ifndef PT_183
    #define PT_183
#endif // PT_183
#ifndef AU_183
    #define AU_183
#endif // AU_183
#ifndef HG_183
    #define HG_183
#endif // HG_183
#ifndef TL_183
    #define TL_183
#endif // TL_183
#ifndef PB_183
    #define PB_183
#endif // PB_183
#ifndef YB_184
    #define YB_184
#endif // YB_184
#ifndef LU_184
    #define LU_184
#endif // LU_184
#ifndef HF_184
    #define HF_184
#endif // HF_184
#ifndef TA_184
    #define TA_184
#endif // TA_184
#ifndef W_184
    #define W_184
#endif // W_184
#ifndef RE_184
    #define RE_184
#endif // RE_184
#ifndef OS_184
    #define OS_184
#endif // OS_184
#ifndef IR_184
    #define IR_184
#endif // IR_184
#ifndef PT_184
    #define PT_184
#endif // PT_184
#ifndef AU_184
    #define AU_184
#endif // AU_184
#ifndef HG_184
    #define HG_184
#endif // HG_184
#ifndef TL_184
    #define TL_184
#endif // TL_184
#ifndef PB_184
    #define PB_184
#endif // PB_184
#ifndef BI_184
    #define BI_184
#endif // BI_184
#ifndef YB_185
    #define YB_185
#endif // YB_185
#ifndef LU_185
    #define LU_185
#endif // LU_185
#ifndef HF_185
    #define HF_185
#endif // HF_185
#ifndef TA_185
    #define TA_185
#endif // TA_185
#ifndef W_185
    #define W_185
#endif // W_185
#ifndef RE_185
    #define RE_185
#endif // RE_185
#ifndef OS_185
    #define OS_185
#endif // OS_185
#ifndef IR_185
    #define IR_185
#endif // IR_185
#ifndef PT_185
    #define PT_185
#endif // PT_185
#ifndef AU_185
    #define AU_185
#endif // AU_185
#ifndef HG_185
    #define HG_185
#endif // HG_185
#ifndef TL_185
    #define TL_185
#endif // TL_185
#ifndef PB_185
    #define PB_185
#endif // PB_185
#ifndef BI_185
    #define BI_185
#endif // BI_185
#ifndef LU_186
    #define LU_186
#endif // LU_186
#ifndef HF_186
    #define HF_186
#endif // HF_186
#ifndef TA_186
    #define TA_186
#endif // TA_186
#ifndef W_186
    #define W_186
#endif // W_186
#ifndef RE_186
    #define RE_186
#endif // RE_186
#ifndef OS_186
    #define OS_186
#endif // OS_186
#ifndef IR_186
    #define IR_186
#endif // IR_186
#ifndef PT_186
    #define PT_186
#endif // PT_186
#ifndef AU_186
    #define AU_186
#endif // AU_186
#ifndef HG_186
    #define HG_186
#endif // HG_186
#ifndef TL_186
    #define TL_186
#endif // TL_186
#ifndef PB_186
    #define PB_186
#endif // PB_186
#ifndef BI_186
    #define BI_186
#endif // BI_186
#ifndef PO_186
    #define PO_186
#endif // PO_186
#ifndef LU_187
    #define LU_187
#endif // LU_187
#ifndef HF_187
    #define HF_187
#endif // HF_187
#ifndef TA_187
    #define TA_187
#endif // TA_187
#ifndef W_187
    #define W_187
#endif // W_187
#ifndef RE_187
    #define RE_187
#endif // RE_187
#ifndef OS_187
    #define OS_187
#endif // OS_187
#ifndef IR_187
    #define IR_187
#endif // IR_187
#ifndef PT_187
    #define PT_187
#endif // PT_187
#ifndef AU_187
    #define AU_187
#endif // AU_187
#ifndef HG_187
    #define HG_187
#endif // HG_187
#ifndef TL_187
    #define TL_187
#endif // TL_187
#ifndef PB_187
    #define PB_187
#endif // PB_187
#ifndef BI_187
    #define BI_187
#endif // BI_187
#ifndef PO_187
    #define PO_187
#endif // PO_187
#ifndef LU_188
    #define LU_188
#endif // LU_188
#ifndef HF_188
    #define HF_188
#endif // HF_188
#ifndef TA_188
    #define TA_188
#endif // TA_188
#ifndef W_188
    #define W_188
#endif // W_188
#ifndef RE_188
    #define RE_188
#endif // RE_188
#ifndef OS_188
    #define OS_188
#endif // OS_188
#ifndef IR_188
    #define IR_188
#endif // IR_188
#ifndef PT_188
    #define PT_188
#endif // PT_188
#ifndef AU_188
    #define AU_188
#endif // AU_188
#ifndef HG_188
    #define HG_188
#endif // HG_188
#ifndef TL_188
    #define TL_188
#endif // TL_188
#ifndef PB_188
    #define PB_188
#endif // PB_188
#ifndef BI_188
    #define BI_188
#endif // BI_188
#ifndef PO_188
    #define PO_188
#endif // PO_188
#ifndef HF_189
    #define HF_189
#endif // HF_189
#ifndef TA_189
    #define TA_189
#endif // TA_189
#ifndef W_189
    #define W_189
#endif // W_189
#ifndef RE_189
    #define RE_189
#endif // RE_189
#ifndef OS_189
    #define OS_189
#endif // OS_189
#ifndef IR_189
    #define IR_189
#endif // IR_189
#ifndef PT_189
    #define PT_189
#endif // PT_189
#ifndef AU_189
    #define AU_189
#endif // AU_189
#ifndef HG_189
    #define HG_189
#endif // HG_189
#ifndef TL_189
    #define TL_189
#endif // TL_189
#ifndef PB_189
    #define PB_189
#endif // PB_189
#ifndef BI_189
    #define BI_189
#endif // BI_189
#ifndef PO_189
    #define PO_189
#endif // PO_189
#ifndef HF_190
    #define HF_190
#endif // HF_190
#ifndef TA_190
    #define TA_190
#endif // TA_190
#ifndef W_190
    #define W_190
#endif // W_190
#ifndef RE_190
    #define RE_190
#endif // RE_190
#ifndef OS_190
    #define OS_190
#endif // OS_190
#ifndef IR_190
    #define IR_190
#endif // IR_190
#ifndef PT_190
    #define PT_190
#endif // PT_190
#ifndef AU_190
    #define AU_190
#endif // AU_190
#ifndef HG_190
    #define HG_190
#endif // HG_190
#ifndef TL_190
    #define TL_190
#endif // TL_190
#ifndef PB_190
    #define PB_190
#endif // PB_190
#ifndef BI_190
    #define BI_190
#endif // BI_190
#ifndef PO_190
    #define PO_190
#endif // PO_190
#ifndef TA_191
    #define TA_191
#endif // TA_191
#ifndef W_191
    #define W_191
#endif // W_191
#ifndef RE_191
    #define RE_191
#endif // RE_191
#ifndef OS_191
    #define OS_191
#endif // OS_191
#ifndef IR_191
    #define IR_191
#endif // IR_191
#ifndef PT_191
    #define PT_191
#endif // PT_191
#ifndef AU_191
    #define AU_191
#endif // AU_191
#ifndef HG_191
    #define HG_191
#endif // HG_191
#ifndef TL_191
    #define TL_191
#endif // TL_191
#ifndef PB_191
    #define PB_191
#endif // PB_191
#ifndef BI_191
    #define BI_191
#endif // BI_191
#ifndef PO_191
    #define PO_191
#endif // PO_191
#ifndef AT_191
    #define AT_191
#endif // AT_191
#ifndef TA_192
    #define TA_192
#endif // TA_192
#ifndef W_192
    #define W_192
#endif // W_192
#ifndef RE_192
    #define RE_192
#endif // RE_192
#ifndef OS_192
    #define OS_192
#endif // OS_192
#ifndef IR_192
    #define IR_192
#endif // IR_192
#ifndef PT_192
    #define PT_192
#endif // PT_192
#ifndef AU_192
    #define AU_192
#endif // AU_192
#ifndef HG_192
    #define HG_192
#endif // HG_192
#ifndef TL_192
    #define TL_192
#endif // TL_192
#ifndef PB_192
    #define PB_192
#endif // PB_192
#ifndef BI_192
    #define BI_192
#endif // BI_192
#ifndef PO_192
    #define PO_192
#endif // PO_192
#ifndef AT_192
    #define AT_192
#endif // AT_192
#ifndef TA_193
    #define TA_193
#endif // TA_193
#ifndef W_193
    #define W_193
#endif // W_193
#ifndef RE_193
    #define RE_193
#endif // RE_193
#ifndef OS_193
    #define OS_193
#endif // OS_193
#ifndef IR_193
    #define IR_193
#endif // IR_193
#ifndef PT_193
    #define PT_193
#endif // PT_193
#ifndef AU_193
    #define AU_193
#endif // AU_193
#ifndef HG_193
    #define HG_193
#endif // HG_193
#ifndef TL_193
    #define TL_193
#endif // TL_193
#ifndef PB_193
    #define PB_193
#endif // PB_193
#ifndef BI_193
    #define BI_193
#endif // BI_193
#ifndef PO_193
    #define PO_193
#endif // PO_193
#ifndef AT_193
    #define AT_193
#endif // AT_193
#ifndef RN_193
    #define RN_193
#endif // RN_193
#ifndef TA_194
    #define TA_194
#endif // TA_194
#ifndef W_194
    #define W_194
#endif // W_194
#ifndef RE_194
    #define RE_194
#endif // RE_194
#ifndef OS_194
    #define OS_194
#endif // OS_194
#ifndef IR_194
    #define IR_194
#endif // IR_194
#ifndef PT_194
    #define PT_194
#endif // PT_194
#ifndef AU_194
    #define AU_194
#endif // AU_194
#ifndef HG_194
    #define HG_194
#endif // HG_194
#ifndef TL_194
    #define TL_194
#endif // TL_194
#ifndef PB_194
    #define PB_194
#endif // PB_194
#ifndef BI_194
    #define BI_194
#endif // BI_194
#ifndef PO_194
    #define PO_194
#endif // PO_194
#ifndef AT_194
    #define AT_194
#endif // AT_194
#ifndef RN_194
    #define RN_194
#endif // RN_194
#ifndef W_195
    #define W_195
#endif // W_195
#ifndef RE_195
    #define RE_195
#endif // RE_195
#ifndef OS_195
    #define OS_195
#endif // OS_195
#ifndef IR_195
    #define IR_195
#endif // IR_195
#ifndef PT_195
    #define PT_195
#endif // PT_195
#ifndef AU_195
    #define AU_195
#endif // AU_195
#ifndef HG_195
    #define HG_195
#endif // HG_195
#ifndef TL_195
    #define TL_195
#endif // TL_195
#ifndef PB_195
    #define PB_195
#endif // PB_195
#ifndef BI_195
    #define BI_195
#endif // BI_195
#ifndef PO_195
    #define PO_195
#endif // PO_195
#ifndef AT_195
    #define AT_195
#endif // AT_195
#ifndef RN_195
    #define RN_195
#endif // RN_195
#ifndef W_196
    #define W_196
#endif // W_196
#ifndef RE_196
    #define RE_196
#endif // RE_196
#ifndef OS_196
    #define OS_196
#endif // OS_196
#ifndef IR_196
    #define IR_196
#endif // IR_196
#ifndef PT_196
    #define PT_196
#endif // PT_196
#ifndef AU_196
    #define AU_196
#endif // AU_196
#ifndef HG_196
    #define HG_196
#endif // HG_196
#ifndef TL_196
    #define TL_196
#endif // TL_196
#ifndef PB_196
    #define PB_196
#endif // PB_196
#ifndef BI_196
    #define BI_196
#endif // BI_196
#ifndef PO_196
    #define PO_196
#endif // PO_196
#ifndef AT_196
    #define AT_196
#endif // AT_196
#ifndef RN_196
    #define RN_196
#endif // RN_196
#ifndef W_197
    #define W_197
#endif // W_197
#ifndef RE_197
    #define RE_197
#endif // RE_197
#ifndef OS_197
    #define OS_197
#endif // OS_197
#ifndef IR_197
    #define IR_197
#endif // IR_197
#ifndef PT_197
    #define PT_197
#endif // PT_197
#ifndef AU_197
    #define AU_197
#endif // AU_197
#ifndef HG_197
    #define HG_197
#endif // HG_197
#ifndef TL_197
    #define TL_197
#endif // TL_197
#ifndef PB_197
    #define PB_197
#endif // PB_197
#ifndef BI_197
    #define BI_197
#endif // BI_197
#ifndef PO_197
    #define PO_197
#endif // PO_197
#ifndef AT_197
    #define AT_197
#endif // AT_197
#ifndef RN_197
    #define RN_197
#endif // RN_197
#ifndef FR_197
    #define FR_197
#endif // FR_197
#ifndef RE_198
    #define RE_198
#endif // RE_198
#ifndef OS_198
    #define OS_198
#endif // OS_198
#ifndef IR_198
    #define IR_198
#endif // IR_198
#ifndef PT_198
    #define PT_198
#endif // PT_198
#ifndef AU_198
    #define AU_198
#endif // AU_198
#ifndef HG_198
    #define HG_198
#endif // HG_198
#ifndef TL_198
    #define TL_198
#endif // TL_198
#ifndef PB_198
    #define PB_198
#endif // PB_198
#ifndef BI_198
    #define BI_198
#endif // BI_198
#ifndef PO_198
    #define PO_198
#endif // PO_198
#ifndef AT_198
    #define AT_198
#endif // AT_198
#ifndef RN_198
    #define RN_198
#endif // RN_198
#ifndef FR_198
    #define FR_198
#endif // FR_198
#ifndef RE_199
    #define RE_199
#endif // RE_199
#ifndef OS_199
    #define OS_199
#endif // OS_199
#ifndef IR_199
    #define IR_199
#endif // IR_199
#ifndef PT_199
    #define PT_199
#endif // PT_199
#ifndef AU_199
    #define AU_199
#endif // AU_199
#ifndef HG_199
    #define HG_199
#endif // HG_199
#ifndef TL_199
    #define TL_199
#endif // TL_199
#ifndef PB_199
    #define PB_199
#endif // PB_199
#ifndef BI_199
    #define BI_199
#endif // BI_199
#ifndef PO_199
    #define PO_199
#endif // PO_199
#ifndef AT_199
    #define AT_199
#endif // AT_199
#ifndef RN_199
    #define RN_199
#endif // RN_199
#ifndef FR_199
    #define FR_199
#endif // FR_199
#ifndef OS_200
    #define OS_200
#endif // OS_200
#ifndef IR_200
    #define IR_200
#endif // IR_200
#ifndef PT_200
    #define PT_200
#endif // PT_200
#ifndef AU_200
    #define AU_200
#endif // AU_200
#ifndef HG_200
    #define HG_200
#endif // HG_200
#ifndef TL_200
    #define TL_200
#endif // TL_200
#ifndef PB_200
    #define PB_200
#endif // PB_200
#ifndef BI_200
    #define BI_200
#endif // BI_200
#ifndef PO_200
    #define PO_200
#endif // PO_200
#ifndef AT_200
    #define AT_200
#endif // AT_200
#ifndef RN_200
    #define RN_200
#endif // RN_200
#ifndef FR_200
    #define FR_200
#endif // FR_200
#ifndef OS_201
    #define OS_201
#endif // OS_201
#ifndef IR_201
    #define IR_201
#endif // IR_201
#ifndef PT_201
    #define PT_201
#endif // PT_201
#ifndef AU_201
    #define AU_201
#endif // AU_201
#ifndef HG_201
    #define HG_201
#endif // HG_201
#ifndef TL_201
    #define TL_201
#endif // TL_201
#ifndef PB_201
    #define PB_201
#endif // PB_201
#ifndef BI_201
    #define BI_201
#endif // BI_201
#ifndef PO_201
    #define PO_201
#endif // PO_201
#ifndef AT_201
    #define AT_201
#endif // AT_201
#ifndef RN_201
    #define RN_201
#endif // RN_201
#ifndef FR_201
    #define FR_201
#endif // FR_201
#ifndef RA_201
    #define RA_201
#endif // RA_201
#ifndef OS_202
    #define OS_202
#endif // OS_202
#ifndef IR_202
    #define IR_202
#endif // IR_202
#ifndef PT_202
    #define PT_202
#endif // PT_202
#ifndef AU_202
    #define AU_202
#endif // AU_202
#ifndef HG_202
    #define HG_202
#endif // HG_202
#ifndef TL_202
    #define TL_202
#endif // TL_202
#ifndef PB_202
    #define PB_202
#endif // PB_202
#ifndef BI_202
    #define BI_202
#endif // BI_202
#ifndef PO_202
    #define PO_202
#endif // PO_202
#ifndef AT_202
    #define AT_202
#endif // AT_202
#ifndef RN_202
    #define RN_202
#endif // RN_202
#ifndef FR_202
    #define FR_202
#endif // FR_202
#ifndef RA_202
    #define RA_202
#endif // RA_202
#ifndef OS_203
    #define OS_203
#endif // OS_203
#ifndef IR_203
    #define IR_203
#endif // IR_203
#ifndef PT_203
    #define PT_203
#endif // PT_203
#ifndef AU_203
    #define AU_203
#endif // AU_203
#ifndef HG_203
    #define HG_203
#endif // HG_203
#ifndef TL_203
    #define TL_203
#endif // TL_203
#ifndef PB_203
    #define PB_203
#endif // PB_203
#ifndef BI_203
    #define BI_203
#endif // BI_203
#ifndef PO_203
    #define PO_203
#endif // PO_203
#ifndef AT_203
    #define AT_203
#endif // AT_203
#ifndef RN_203
    #define RN_203
#endif // RN_203
#ifndef FR_203
    #define FR_203
#endif // FR_203
#ifndef RA_203
    #define RA_203
#endif // RA_203
#ifndef IR_204
    #define IR_204
#endif // IR_204
#ifndef PT_204
    #define PT_204
#endif // PT_204
#ifndef AU_204
    #define AU_204
#endif // AU_204
#ifndef HG_204
    #define HG_204
#endif // HG_204
#ifndef TL_204
    #define TL_204
#endif // TL_204
#ifndef PB_204
    #define PB_204
#endif // PB_204
#ifndef BI_204
    #define BI_204
#endif // BI_204
#ifndef PO_204
    #define PO_204
#endif // PO_204
#ifndef AT_204
    #define AT_204
#endif // AT_204
#ifndef RN_204
    #define RN_204
#endif // RN_204
#ifndef FR_204
    #define FR_204
#endif // FR_204
#ifndef RA_204
    #define RA_204
#endif // RA_204
#ifndef IR_205
    #define IR_205
#endif // IR_205
#ifndef PT_205
    #define PT_205
#endif // PT_205
#ifndef AU_205
    #define AU_205
#endif // AU_205
#ifndef HG_205
    #define HG_205
#endif // HG_205
#ifndef TL_205
    #define TL_205
#endif // TL_205
#ifndef PB_205
    #define PB_205
#endif // PB_205
#ifndef BI_205
    #define BI_205
#endif // BI_205
#ifndef PO_205
    #define PO_205
#endif // PO_205
#ifndef AT_205
    #define AT_205
#endif // AT_205
#ifndef RN_205
    #define RN_205
#endif // RN_205
#ifndef FR_205
    #define FR_205
#endif // FR_205
#ifndef RA_205
    #define RA_205
#endif // RA_205
#ifndef AC_205
    #define AC_205
#endif // AC_205
#ifndef PT_206
    #define PT_206
#endif // PT_206
#ifndef AU_206
    #define AU_206
#endif // AU_206
#ifndef HG_206
    #define HG_206
#endif // HG_206
#ifndef TL_206
    #define TL_206
#endif // TL_206
#ifndef PB_206
    #define PB_206
#endif // PB_206
#ifndef BI_206
    #define BI_206
#endif // BI_206
#ifndef PO_206
    #define PO_206
#endif // PO_206
#ifndef AT_206
    #define AT_206
#endif // AT_206
#ifndef RN_206
    #define RN_206
#endif // RN_206
#ifndef FR_206
    #define FR_206
#endif // FR_206
#ifndef RA_206
    #define RA_206
#endif // RA_206
#ifndef AC_206
    #define AC_206
#endif // AC_206
#ifndef PT_207
    #define PT_207
#endif // PT_207
#ifndef AU_207
    #define AU_207
#endif // AU_207
#ifndef HG_207
    #define HG_207
#endif // HG_207
#ifndef TL_207
    #define TL_207
#endif // TL_207
#ifndef PB_207
    #define PB_207
#endif // PB_207
#ifndef BI_207
    #define BI_207
#endif // BI_207
#ifndef PO_207
    #define PO_207
#endif // PO_207
#ifndef AT_207
    #define AT_207
#endif // AT_207
#ifndef RN_207
    #define RN_207
#endif // RN_207
#ifndef FR_207
    #define FR_207
#endif // FR_207
#ifndef RA_207
    #define RA_207
#endif // RA_207
#ifndef AC_207
    #define AC_207
#endif // AC_207
#ifndef PT_208
    #define PT_208
#endif // PT_208
#ifndef AU_208
    #define AU_208
#endif // AU_208
#ifndef HG_208
    #define HG_208
#endif // HG_208
#ifndef TL_208
    #define TL_208
#endif // TL_208
#ifndef PB_208
    #define PB_208
#endif // PB_208
#ifndef BI_208
    #define BI_208
#endif // BI_208
#ifndef PO_208
    #define PO_208
#endif // PO_208
#ifndef AT_208
    #define AT_208
#endif // AT_208
#ifndef RN_208
    #define RN_208
#endif // RN_208
#ifndef FR_208
    #define FR_208
#endif // FR_208
#ifndef RA_208
    #define RA_208
#endif // RA_208
#ifndef AC_208
    #define AC_208
#endif // AC_208
#ifndef TH_208
    #define TH_208
#endif // TH_208
#ifndef AU_209
    #define AU_209
#endif // AU_209
#ifndef HG_209
    #define HG_209
#endif // HG_209
#ifndef TL_209
    #define TL_209
#endif // TL_209
#ifndef PB_209
    #define PB_209
#endif // PB_209
#ifndef BI_209
    #define BI_209
#endif // BI_209
#ifndef PO_209
    #define PO_209
#endif // PO_209
#ifndef AT_209
    #define AT_209
#endif // AT_209
#ifndef RN_209
    #define RN_209
#endif // RN_209
#ifndef FR_209
    #define FR_209
#endif // FR_209
#ifndef RA_209
    #define RA_209
#endif // RA_209
#ifndef AC_209
    #define AC_209
#endif // AC_209
#ifndef TH_209
    #define TH_209
#endif // TH_209
#ifndef AU_210
    #define AU_210
#endif // AU_210
#ifndef HG_210
    #define HG_210
#endif // HG_210
#ifndef TL_210
    #define TL_210
#endif // TL_210
#ifndef PB_210
    #define PB_210
#endif // PB_210
#ifndef BI_210
    #define BI_210
#endif // BI_210
#ifndef PO_210
    #define PO_210
#endif // PO_210
#ifndef AT_210
    #define AT_210
#endif // AT_210
#ifndef RN_210
    #define RN_210
#endif // RN_210
#ifndef FR_210
    #define FR_210
#endif // FR_210
#ifndef RA_210
    #define RA_210
#endif // RA_210
#ifndef AC_210
    #define AC_210
#endif // AC_210
#ifndef TH_210
    #define TH_210
#endif // TH_210
#ifndef HG_211
    #define HG_211
#endif // HG_211
#ifndef TL_211
    #define TL_211
#endif // TL_211
#ifndef PB_211
    #define PB_211
#endif // PB_211
#ifndef BI_211
    #define BI_211
#endif // BI_211
#ifndef PO_211
    #define PO_211
#endif // PO_211
#ifndef AT_211
    #define AT_211
#endif // AT_211
#ifndef RN_211
    #define RN_211
#endif // RN_211
#ifndef FR_211
    #define FR_211
#endif // FR_211
#ifndef RA_211
    #define RA_211
#endif // RA_211
#ifndef AC_211
    #define AC_211
#endif // AC_211
#ifndef TH_211
    #define TH_211
#endif // TH_211
#ifndef PA_211
    #define PA_211
#endif // PA_211
#ifndef HG_212
    #define HG_212
#endif // HG_212
#ifndef TL_212
    #define TL_212
#endif // TL_212
#ifndef PB_212
    #define PB_212
#endif // PB_212
#ifndef BI_212
    #define BI_212
#endif // BI_212
#ifndef PO_212
    #define PO_212
#endif // PO_212
#ifndef AT_212
    #define AT_212
#endif // AT_212
#ifndef RN_212
    #define RN_212
#endif // RN_212
#ifndef FR_212
    #define FR_212
#endif // FR_212
#ifndef RA_212
    #define RA_212
#endif // RA_212
#ifndef AC_212
    #define AC_212
#endif // AC_212
#ifndef TH_212
    #define TH_212
#endif // TH_212
#ifndef PA_212
    #define PA_212
#endif // PA_212
#ifndef HG_213
    #define HG_213
#endif // HG_213
#ifndef TL_213
    #define TL_213
#endif // TL_213
#ifndef PB_213
    #define PB_213
#endif // PB_213
#ifndef BI_213
    #define BI_213
#endif // BI_213
#ifndef PO_213
    #define PO_213
#endif // PO_213
#ifndef AT_213
    #define AT_213
#endif // AT_213
#ifndef RN_213
    #define RN_213
#endif // RN_213
#ifndef FR_213
    #define FR_213
#endif // FR_213
#ifndef RA_213
    #define RA_213
#endif // RA_213
#ifndef AC_213
    #define AC_213
#endif // AC_213
#ifndef TH_213
    #define TH_213
#endif // TH_213
#ifndef PA_213
    #define PA_213
#endif // PA_213
#ifndef HG_214
    #define HG_214
#endif // HG_214
#ifndef TL_214
    #define TL_214
#endif // TL_214
#ifndef PB_214
    #define PB_214
#endif // PB_214
#ifndef BI_214
    #define BI_214
#endif // BI_214
#ifndef PO_214
    #define PO_214
#endif // PO_214
#ifndef AT_214
    #define AT_214
#endif // AT_214
#ifndef RN_214
    #define RN_214
#endif // RN_214
#ifndef FR_214
    #define FR_214
#endif // FR_214
#ifndef RA_214
    #define RA_214
#endif // RA_214
#ifndef AC_214
    #define AC_214
#endif // AC_214
#ifndef TH_214
    #define TH_214
#endif // TH_214
#ifndef PA_214
    #define PA_214
#endif // PA_214
#ifndef HG_215
    #define HG_215
#endif // HG_215
#ifndef TL_215
    #define TL_215
#endif // TL_215
#ifndef PB_215
    #define PB_215
#endif // PB_215
#ifndef BI_215
    #define BI_215
#endif // BI_215
#ifndef PO_215
    #define PO_215
#endif // PO_215
#ifndef AT_215
    #define AT_215
#endif // AT_215
#ifndef RN_215
    #define RN_215
#endif // RN_215
#ifndef FR_215
    #define FR_215
#endif // FR_215
#ifndef RA_215
    #define RA_215
#endif // RA_215
#ifndef AC_215
    #define AC_215
#endif // AC_215
#ifndef TH_215
    #define TH_215
#endif // TH_215
#ifndef PA_215
    #define PA_215
#endif // PA_215
#ifndef U_215
    #define U_215
#endif // U_215
#ifndef HG_216
    #define HG_216
#endif // HG_216
#ifndef TL_216
    #define TL_216
#endif // TL_216
#ifndef PB_216
    #define PB_216
#endif // PB_216
#ifndef BI_216
    #define BI_216
#endif // BI_216
#ifndef PO_216
    #define PO_216
#endif // PO_216
#ifndef AT_216
    #define AT_216
#endif // AT_216
#ifndef RN_216
    #define RN_216
#endif // RN_216
#ifndef FR_216
    #define FR_216
#endif // FR_216
#ifndef RA_216
    #define RA_216
#endif // RA_216
#ifndef AC_216
    #define AC_216
#endif // AC_216
#ifndef TH_216
    #define TH_216
#endif // TH_216
#ifndef PA_216
    #define PA_216
#endif // PA_216
#ifndef U_216
    #define U_216
#endif // U_216
#ifndef TL_217
    #define TL_217
#endif // TL_217
#ifndef PB_217
    #define PB_217
#endif // PB_217
#ifndef BI_217
    #define BI_217
#endif // BI_217
#ifndef PO_217
    #define PO_217
#endif // PO_217
#ifndef AT_217
    #define AT_217
#endif // AT_217
#ifndef RN_217
    #define RN_217
#endif // RN_217
#ifndef FR_217
    #define FR_217
#endif // FR_217
#ifndef RA_217
    #define RA_217
#endif // RA_217
#ifndef AC_217
    #define AC_217
#endif // AC_217
#ifndef TH_217
    #define TH_217
#endif // TH_217
#ifndef PA_217
    #define PA_217
#endif // PA_217
#ifndef U_217
    #define U_217
#endif // U_217
#ifndef TL_218
    #define TL_218
#endif // TL_218
#ifndef PB_218
    #define PB_218
#endif // PB_218
#ifndef BI_218
    #define BI_218
#endif // BI_218
#ifndef PO_218
    #define PO_218
#endif // PO_218
#ifndef AT_218
    #define AT_218
#endif // AT_218
#ifndef RN_218
    #define RN_218
#endif // RN_218
#ifndef FR_218
    #define FR_218
#endif // FR_218
#ifndef RA_218
    #define RA_218
#endif // RA_218
#ifndef AC_218
    #define AC_218
#endif // AC_218
#ifndef TH_218
    #define TH_218
#endif // TH_218
#ifndef PA_218
    #define PA_218
#endif // PA_218
#ifndef U_218
    #define U_218
#endif // U_218
#ifndef PB_219
    #define PB_219
#endif // PB_219
#ifndef BI_219
    #define BI_219
#endif // BI_219
#ifndef PO_219
    #define PO_219
#endif // PO_219
#ifndef AT_219
    #define AT_219
#endif // AT_219
#ifndef RN_219
    #define RN_219
#endif // RN_219
#ifndef FR_219
    #define FR_219
#endif // FR_219
#ifndef RA_219
    #define RA_219
#endif // RA_219
#ifndef AC_219
    #define AC_219
#endif // AC_219
#ifndef TH_219
    #define TH_219
#endif // TH_219
#ifndef PA_219
    #define PA_219
#endif // PA_219
#ifndef U_219
    #define U_219
#endif // U_219
#ifndef NP_219
    #define NP_219
#endif // NP_219
#ifndef PB_220
    #define PB_220
#endif // PB_220
#ifndef BI_220
    #define BI_220
#endif // BI_220
#ifndef PO_220
    #define PO_220
#endif // PO_220
#ifndef AT_220
    #define AT_220
#endif // AT_220
#ifndef RN_220
    #define RN_220
#endif // RN_220
#ifndef FR_220
    #define FR_220
#endif // FR_220
#ifndef RA_220
    #define RA_220
#endif // RA_220
#ifndef AC_220
    #define AC_220
#endif // AC_220
#ifndef TH_220
    #define TH_220
#endif // TH_220
#ifndef PA_220
    #define PA_220
#endif // PA_220
#ifndef U_220
    #define U_220
#endif // U_220
#ifndef NP_220
    #define NP_220
#endif // NP_220
#ifndef BI_221
    #define BI_221
#endif // BI_221
#ifndef PO_221
    #define PO_221
#endif // PO_221
#ifndef AT_221
    #define AT_221
#endif // AT_221
#ifndef RN_221
    #define RN_221
#endif // RN_221
#ifndef FR_221
    #define FR_221
#endif // FR_221
#ifndef RA_221
    #define RA_221
#endif // RA_221
#ifndef AC_221
    #define AC_221
#endif // AC_221
#ifndef TH_221
    #define TH_221
#endif // TH_221
#ifndef PA_221
    #define PA_221
#endif // PA_221
#ifndef U_221
    #define U_221
#endif // U_221
#ifndef NP_221
    #define NP_221
#endif // NP_221
#ifndef PU_221
    #define PU_221
#endif // PU_221
#ifndef BI_222
    #define BI_222
#endif // BI_222
#ifndef PO_222
    #define PO_222
#endif // PO_222
#ifndef AT_222
    #define AT_222
#endif // AT_222
#ifndef RN_222
    #define RN_222
#endif // RN_222
#ifndef FR_222
    #define FR_222
#endif // FR_222
#ifndef RA_222
    #define RA_222
#endif // RA_222
#ifndef AC_222
    #define AC_222
#endif // AC_222
#ifndef TH_222
    #define TH_222
#endif // TH_222
#ifndef PA_222
    #define PA_222
#endif // PA_222
#ifndef U_222
    #define U_222
#endif // U_222
#ifndef NP_222
    #define NP_222
#endif // NP_222
#ifndef PU_222
    #define PU_222
#endif // PU_222
#ifndef BI_223
    #define BI_223
#endif // BI_223
#ifndef PO_223
    #define PO_223
#endif // PO_223
#ifndef AT_223
    #define AT_223
#endif // AT_223
#ifndef RN_223
    #define RN_223
#endif // RN_223
#ifndef FR_223
    #define FR_223
#endif // FR_223
#ifndef RA_223
    #define RA_223
#endif // RA_223
#ifndef AC_223
    #define AC_223
#endif // AC_223
#ifndef TH_223
    #define TH_223
#endif // TH_223
#ifndef PA_223
    #define PA_223
#endif // PA_223
#ifndef U_223
    #define U_223
#endif // U_223
#ifndef NP_223
    #define NP_223
#endif // NP_223
#ifndef PU_223
    #define PU_223
#endif // PU_223
#ifndef AM_223
    #define AM_223
#endif // AM_223
#ifndef BI_224
    #define BI_224
#endif // BI_224
#ifndef PO_224
    #define PO_224
#endif // PO_224
#ifndef AT_224
    #define AT_224
#endif // AT_224
#ifndef RN_224
    #define RN_224
#endif // RN_224
#ifndef FR_224
    #define FR_224
#endif // FR_224
#ifndef RA_224
    #define RA_224
#endif // RA_224
#ifndef AC_224
    #define AC_224
#endif // AC_224
#ifndef TH_224
    #define TH_224
#endif // TH_224
#ifndef PA_224
    #define PA_224
#endif // PA_224
#ifndef U_224
    #define U_224
#endif // U_224
#ifndef NP_224
    #define NP_224
#endif // NP_224
#ifndef PU_224
    #define PU_224
#endif // PU_224
#ifndef AM_224
    #define AM_224
#endif // AM_224
#ifndef PO_225
    #define PO_225
#endif // PO_225
#ifndef AT_225
    #define AT_225
#endif // AT_225
#ifndef RN_225
    #define RN_225
#endif // RN_225
#ifndef FR_225
    #define FR_225
#endif // FR_225
#ifndef RA_225
    #define RA_225
#endif // RA_225
#ifndef AC_225
    #define AC_225
#endif // AC_225
#ifndef TH_225
    #define TH_225
#endif // TH_225
#ifndef PA_225
    #define PA_225
#endif // PA_225
#ifndef U_225
    #define U_225
#endif // U_225
#ifndef NP_225
    #define NP_225
#endif // NP_225
#ifndef PU_225
    #define PU_225
#endif // PU_225
#ifndef AM_225
    #define AM_225
#endif // AM_225
#ifndef PO_226
    #define PO_226
#endif // PO_226
#ifndef AT_226
    #define AT_226
#endif // AT_226
#ifndef RN_226
    #define RN_226
#endif // RN_226
#ifndef FR_226
    #define FR_226
#endif // FR_226
#ifndef RA_226
    #define RA_226
#endif // RA_226
#ifndef AC_226
    #define AC_226
#endif // AC_226
#ifndef TH_226
    #define TH_226
#endif // TH_226
#ifndef PA_226
    #define PA_226
#endif // PA_226
#ifndef U_226
    #define U_226
#endif // U_226
#ifndef NP_226
    #define NP_226
#endif // NP_226
#ifndef PU_226
    #define PU_226
#endif // PU_226
#ifndef AM_226
    #define AM_226
#endif // AM_226
#ifndef PO_227
    #define PO_227
#endif // PO_227
#ifndef AT_227
    #define AT_227
#endif // AT_227
#ifndef RN_227
    #define RN_227
#endif // RN_227
#ifndef FR_227
    #define FR_227
#endif // FR_227
#ifndef RA_227
    #define RA_227
#endif // RA_227
#ifndef AC_227
    #define AC_227
#endif // AC_227
#ifndef TH_227
    #define TH_227
#endif // TH_227
#ifndef PA_227
    #define PA_227
#endif // PA_227
#ifndef U_227
    #define U_227
#endif // U_227
#ifndef NP_227
    #define NP_227
#endif // NP_227
#ifndef PU_227
    #define PU_227
#endif // PU_227
#ifndef AM_227
    #define AM_227
#endif // AM_227
#ifndef AT_228
    #define AT_228
#endif // AT_228
#ifndef RN_228
    #define RN_228
#endif // RN_228
#ifndef FR_228
    #define FR_228
#endif // FR_228
#ifndef RA_228
    #define RA_228
#endif // RA_228
#ifndef AC_228
    #define AC_228
#endif // AC_228
#ifndef TH_228
    #define TH_228
#endif // TH_228
#ifndef PA_228
    #define PA_228
#endif // PA_228
#ifndef U_228
    #define U_228
#endif // U_228
#ifndef NP_228
    #define NP_228
#endif // NP_228
#ifndef PU_228
    #define PU_228
#endif // PU_228
#ifndef AM_228
    #define AM_228
#endif // AM_228
#ifndef AT_229
    #define AT_229
#endif // AT_229
#ifndef RN_229
    #define RN_229
#endif // RN_229
#ifndef FR_229
    #define FR_229
#endif // FR_229
#ifndef RA_229
    #define RA_229
#endif // RA_229
#ifndef AC_229
    #define AC_229
#endif // AC_229
#ifndef TH_229
    #define TH_229
#endif // TH_229
#ifndef PA_229
    #define PA_229
#endif // PA_229
#ifndef U_229
    #define U_229
#endif // U_229
#ifndef NP_229
    #define NP_229
#endif // NP_229
#ifndef PU_229
    #define PU_229
#endif // PU_229
#ifndef AM_229
    #define AM_229
#endif // AM_229
#ifndef RN_230
    #define RN_230
#endif // RN_230
#ifndef FR_230
    #define FR_230
#endif // FR_230
#ifndef RA_230
    #define RA_230
#endif // RA_230
#ifndef AC_230
    #define AC_230
#endif // AC_230
#ifndef TH_230
    #define TH_230
#endif // TH_230
#ifndef PA_230
    #define PA_230
#endif // PA_230
#ifndef U_230
    #define U_230
#endif // U_230
#ifndef NP_230
    #define NP_230
#endif // NP_230
#ifndef PU_230
    #define PU_230
#endif // PU_230
#ifndef AM_230
    #define AM_230
#endif // AM_230
#ifndef RN_231
    #define RN_231
#endif // RN_231
#ifndef FR_231
    #define FR_231
#endif // FR_231
#ifndef RA_231
    #define RA_231
#endif // RA_231
#ifndef AC_231
    #define AC_231
#endif // AC_231
#ifndef TH_231
    #define TH_231
#endif // TH_231
#ifndef PA_231
    #define PA_231
#endif // PA_231
#ifndef U_231
    #define U_231
#endif // U_231
#ifndef NP_231
    #define NP_231
#endif // NP_231
#ifndef PU_231
    #define PU_231
#endif // PU_231
#ifndef AM_231
    #define AM_231
#endif // AM_231
#ifndef CM_231
    #define CM_231
#endif // CM_231
#ifndef FR_232
    #define FR_232
#endif // FR_232
#ifndef RA_232
    #define RA_232
#endif // RA_232
#ifndef AC_232
    #define AC_232
#endif // AC_232
#ifndef TH_232
    #define TH_232
#endif // TH_232
#ifndef PA_232
    #define PA_232
#endif // PA_232
#ifndef U_232
    #define U_232
#endif // U_232
#ifndef NP_232
    #define NP_232
#endif // NP_232
#ifndef PU_232
    #define PU_232
#endif // PU_232
#ifndef AM_232
    #define AM_232
#endif // AM_232
#ifndef CM_232
    #define CM_232
#endif // CM_232
#ifndef FR_233
    #define FR_233
#endif // FR_233
#ifndef RA_233
    #define RA_233
#endif // RA_233
#ifndef AC_233
    #define AC_233
#endif // AC_233
#ifndef TH_233
    #define TH_233
#endif // TH_233
#ifndef PA_233
    #define PA_233
#endif // PA_233
#ifndef U_233
    #define U_233
#endif // U_233
#ifndef NP_233
    #define NP_233
#endif // NP_233
#ifndef PU_233
    #define PU_233
#endif // PU_233
#ifndef AM_233
    #define AM_233
#endif // AM_233
#ifndef CM_233
    #define CM_233
#endif // CM_233
#ifndef BK_233
    #define BK_233
#endif // BK_233
#ifndef RA_234
    #define RA_234
#endif // RA_234
#ifndef AC_234
    #define AC_234
#endif // AC_234
#ifndef TH_234
    #define TH_234
#endif // TH_234
#ifndef PA_234
    #define PA_234
#endif // PA_234
#ifndef U_234
    #define U_234
#endif // U_234
#ifndef NP_234
    #define NP_234
#endif // NP_234
#ifndef PU_234
    #define PU_234
#endif // PU_234
#ifndef AM_234
    #define AM_234
#endif // AM_234
#ifndef CM_234
    #define CM_234
#endif // CM_234
#ifndef BK_234
    #define BK_234
#endif // BK_234
#ifndef RA_235
    #define RA_235
#endif // RA_235
#ifndef AC_235
    #define AC_235
#endif // AC_235
#ifndef TH_235
    #define TH_235
#endif // TH_235
#ifndef PA_235
    #define PA_235
#endif // PA_235
#ifndef U_235
    #define U_235
#endif // U_235
#ifndef NP_235
    #define NP_235
#endif // NP_235
#ifndef PU_235
    #define PU_235
#endif // PU_235
#ifndef AM_235
    #define AM_235
#endif // AM_235
#ifndef CM_235
    #define CM_235
#endif // CM_235
#ifndef BK_235
    #define BK_235
#endif // BK_235
#ifndef AC_236
    #define AC_236
#endif // AC_236
#ifndef TH_236
    #define TH_236
#endif // TH_236
#ifndef PA_236
    #define PA_236
#endif // PA_236
#ifndef U_236
    #define U_236
#endif // U_236
#ifndef NP_236
    #define NP_236
#endif // NP_236
#ifndef PU_236
    #define PU_236
#endif // PU_236
#ifndef AM_236
    #define AM_236
#endif // AM_236
#ifndef CM_236
    #define CM_236
#endif // CM_236
#ifndef BK_236
    #define BK_236
#endif // BK_236
#ifndef AC_237
    #define AC_237
#endif // AC_237
#ifndef TH_237
    #define TH_237
#endif // TH_237
#ifndef PA_237
    #define PA_237
#endif // PA_237
#ifndef U_237
    #define U_237
#endif // U_237
#ifndef NP_237
    #define NP_237
#endif // NP_237
#ifndef PU_237
    #define PU_237
#endif // PU_237
#ifndef AM_237
    #define AM_237
#endif // AM_237
#ifndef CM_237
    #define CM_237
#endif // CM_237
#ifndef BK_237
    #define BK_237
#endif // BK_237
#ifndef CF_237
    #define CF_237
#endif // CF_237
#ifndef TH_238
    #define TH_238
#endif // TH_238
#ifndef PA_238
    #define PA_238
#endif // PA_238
#ifndef U_238
    #define U_238
#endif // U_238
#ifndef NP_238
    #define NP_238
#endif // NP_238
#ifndef PU_238
    #define PU_238
#endif // PU_238
#ifndef AM_238
    #define AM_238
#endif // AM_238
#ifndef CM_238
    #define CM_238
#endif // CM_238
#ifndef BK_238
    #define BK_238
#endif // BK_238
#ifndef CF_238
    #define CF_238
#endif // CF_238
#ifndef TH_239
    #define TH_239
#endif // TH_239
#ifndef PA_239
    #define PA_239
#endif // PA_239
#ifndef U_239
    #define U_239
#endif // U_239
#ifndef NP_239
    #define NP_239
#endif // NP_239
#ifndef PU_239
    #define PU_239
#endif // PU_239
#ifndef AM_239
    #define AM_239
#endif // AM_239
#ifndef CM_239
    #define CM_239
#endif // CM_239
#ifndef BK_239
    #define BK_239
#endif // BK_239
#ifndef CF_239
    #define CF_239
#endif // CF_239
#ifndef ES_239
    #define ES_239
#endif // ES_239
#ifndef PA_240
    #define PA_240
#endif // PA_240
#ifndef U_240
    #define U_240
#endif // U_240
#ifndef NP_240
    #define NP_240
#endif // NP_240
#ifndef PU_240
    #define PU_240
#endif // PU_240
#ifndef AM_240
    #define AM_240
#endif // AM_240
#ifndef CM_240
    #define CM_240
#endif // CM_240
#ifndef BK_240
    #define BK_240
#endif // BK_240
#ifndef CF_240
    #define CF_240
#endif // CF_240
#ifndef ES_240
    #define ES_240
#endif // ES_240
#ifndef PA_241
    #define PA_241
#endif // PA_241
#ifndef U_241
    #define U_241
#endif // U_241
#ifndef NP_241
    #define NP_241
#endif // NP_241
#ifndef PU_241
    #define PU_241
#endif // PU_241
#ifndef AM_241
    #define AM_241
#endif // AM_241
#ifndef CM_241
    #define CM_241
#endif // CM_241
#ifndef BK_241
    #define BK_241
#endif // BK_241
#ifndef CF_241
    #define CF_241
#endif // CF_241
#ifndef ES_241
    #define ES_241
#endif // ES_241
#ifndef FM_241
    #define FM_241
#endif // FM_241
#ifndef U_242
    #define U_242
#endif // U_242
#ifndef NP_242
    #define NP_242
#endif // NP_242
#ifndef PU_242
    #define PU_242
#endif // PU_242
#ifndef AM_242
    #define AM_242
#endif // AM_242
#ifndef CM_242
    #define CM_242
#endif // CM_242
#ifndef BK_242
    #define BK_242
#endif // BK_242
#ifndef CF_242
    #define CF_242
#endif // CF_242
#ifndef ES_242
    #define ES_242
#endif // ES_242
#ifndef FM_242
    #define FM_242
#endif // FM_242
#ifndef U_243
    #define U_243
#endif // U_243
#ifndef NP_243
    #define NP_243
#endif // NP_243
#ifndef PU_243
    #define PU_243
#endif // PU_243
#ifndef AM_243
    #define AM_243
#endif // AM_243
#ifndef CM_243
    #define CM_243
#endif // CM_243
#ifndef BK_243
    #define BK_243
#endif // BK_243
#ifndef CF_243
    #define CF_243
#endif // CF_243
#ifndef ES_243
    #define ES_243
#endif // ES_243
#ifndef FM_243
    #define FM_243
#endif // FM_243
#ifndef NP_244
    #define NP_244
#endif // NP_244
#ifndef PU_244
    #define PU_244
#endif // PU_244
#ifndef AM_244
    #define AM_244
#endif // AM_244
#ifndef CM_244
    #define CM_244
#endif // CM_244
#ifndef BK_244
    #define BK_244
#endif // BK_244
#ifndef CF_244
    #define CF_244
#endif // CF_244
#ifndef ES_244
    #define ES_244
#endif // ES_244
#ifndef FM_244
    #define FM_244
#endif // FM_244
#ifndef MD_244
    #define MD_244
#endif // MD_244
#ifndef NP_245
    #define NP_245
#endif // NP_245
#ifndef PU_245
    #define PU_245
#endif // PU_245
#ifndef AM_245
    #define AM_245
#endif // AM_245
#ifndef CM_245
    #define CM_245
#endif // CM_245
#ifndef BK_245
    #define BK_245
#endif // BK_245
#ifndef CF_245
    #define CF_245
#endif // CF_245
#ifndef ES_245
    #define ES_245
#endif // ES_245
#ifndef FM_245
    #define FM_245
#endif // FM_245
#ifndef MD_245
    #define MD_245
#endif // MD_245
#ifndef PU_246
    #define PU_246
#endif // PU_246
#ifndef AM_246
    #define AM_246
#endif // AM_246
#ifndef CM_246
    #define CM_246
#endif // CM_246
#ifndef BK_246
    #define BK_246
#endif // BK_246
#ifndef CF_246
    #define CF_246
#endif // CF_246
#ifndef ES_246
    #define ES_246
#endif // ES_246
#ifndef FM_246
    #define FM_246
#endif // FM_246
#ifndef MD_246
    #define MD_246
#endif // MD_246
#ifndef PU_247
    #define PU_247
#endif // PU_247
#ifndef AM_247
    #define AM_247
#endif // AM_247
#ifndef CM_247
    #define CM_247
#endif // CM_247
#ifndef BK_247
    #define BK_247
#endif // BK_247
#ifndef CF_247
    #define CF_247
#endif // CF_247
#ifndef ES_247
    #define ES_247
#endif // ES_247
#ifndef FM_247
    #define FM_247
#endif // FM_247
#ifndef MD_247
    #define MD_247
#endif // MD_247
#ifndef AM_248
    #define AM_248
#endif // AM_248
#ifndef CM_248
    #define CM_248
#endif // CM_248
#ifndef BK_248
    #define BK_248
#endif // BK_248
#ifndef CF_248
    #define CF_248
#endif // CF_248
#ifndef ES_248
    #define ES_248
#endif // ES_248
#ifndef FM_248
    #define FM_248
#endif // FM_248
#ifndef MD_248
    #define MD_248
#endif // MD_248
#ifndef NO_248
    #define NO_248
#endif // NO_248
#ifndef AM_249
    #define AM_249
#endif // AM_249
#ifndef CM_249
    #define CM_249
#endif // CM_249
#ifndef BK_249
    #define BK_249
#endif // BK_249
#ifndef CF_249
    #define CF_249
#endif // CF_249
#ifndef ES_249
    #define ES_249
#endif // ES_249
#ifndef FM_249
    #define FM_249
#endif // FM_249
#ifndef MD_249
    #define MD_249
#endif // MD_249
#ifndef NO_249
    #define NO_249
#endif // NO_249
#ifndef CM_250
    #define CM_250
#endif // CM_250
#ifndef BK_250
    #define BK_250
#endif // BK_250
#ifndef CF_250
    #define CF_250
#endif // CF_250
#ifndef ES_250
    #define ES_250
#endif // ES_250
#ifndef FM_250
    #define FM_250
#endif // FM_250
#ifndef MD_250
    #define MD_250
#endif // MD_250
#ifndef NO_250
    #define NO_250
#endif // NO_250
#ifndef CM_251
    #define CM_251
#endif // CM_251
#ifndef BK_251
    #define BK_251
#endif // BK_251
#ifndef CF_251
    #define CF_251
#endif // CF_251
#ifndef ES_251
    #define ES_251
#endif // ES_251
#ifndef FM_251
    #define FM_251
#endif // FM_251
#ifndef MD_251
    #define MD_251
#endif // MD_251
#ifndef NO_251
    #define NO_251
#endif // NO_251
#ifndef LR_251
    #define LR_251
#endif // LR_251
#ifndef CM_252
    #define CM_252
#endif // CM_252
#ifndef BK_252
    #define BK_252
#endif // BK_252
#ifndef CF_252
    #define CF_252
#endif // CF_252
#ifndef ES_252
    #define ES_252
#endif // ES_252
#ifndef FM_252
    #define FM_252
#endif // FM_252
#ifndef MD_252
    #define MD_252
#endif // MD_252
#ifndef NO_252
    #define NO_252
#endif // NO_252
#ifndef LR_252
    #define LR_252
#endif // LR_252
#ifndef BK_253
    #define BK_253
#endif // BK_253
#ifndef CF_253
    #define CF_253
#endif // CF_253
#ifndef ES_253
    #define ES_253
#endif // ES_253
#ifndef FM_253
    #define FM_253
#endif // FM_253
#ifndef MD_253
    #define MD_253
#endif // MD_253
#ifndef NO_253
    #define NO_253
#endif // NO_253
#ifndef LR_253
    #define LR_253
#endif // LR_253
#ifndef RF_253
    #define RF_253
#endif // RF_253
#ifndef BK_254
    #define BK_254
#endif // BK_254
#ifndef CF_254
    #define CF_254
#endif // CF_254
#ifndef ES_254
    #define ES_254
#endif // ES_254
#ifndef FM_254
    #define FM_254
#endif // FM_254
#ifndef MD_254
    #define MD_254
#endif // MD_254
#ifndef NO_254
    #define NO_254
#endif // NO_254
#ifndef LR_254
    #define LR_254
#endif // LR_254
#ifndef RF_254
    #define RF_254
#endif // RF_254
#ifndef CF_255
    #define CF_255
#endif // CF_255
#ifndef ES_255
    #define ES_255
#endif // ES_255
#ifndef FM_255
    #define FM_255
#endif // FM_255
#ifndef MD_255
    #define MD_255
#endif // MD_255
#ifndef NO_255
    #define NO_255
#endif // NO_255
#ifndef LR_255
    #define LR_255
#endif // LR_255
#ifndef RF_255
    #define RF_255
#endif // RF_255
#ifndef DB_255
    #define DB_255
#endif // DB_255
#ifndef CF_256
    #define CF_256
#endif // CF_256
#ifndef ES_256
    #define ES_256
#endif // ES_256
#ifndef FM_256
    #define FM_256
#endif // FM_256
#ifndef MD_256
    #define MD_256
#endif // MD_256
#ifndef NO_256
    #define NO_256
#endif // NO_256
#ifndef LR_256
    #define LR_256
#endif // LR_256
#ifndef RF_256
    #define RF_256
#endif // RF_256
#ifndef DB_256
    #define DB_256
#endif // DB_256
#ifndef ES_257
    #define ES_257
#endif // ES_257
#ifndef FM_257
    #define FM_257
#endif // FM_257
#ifndef MD_257
    #define MD_257
#endif // MD_257
#ifndef NO_257
    #define NO_257
#endif // NO_257
#ifndef LR_257
    #define LR_257
#endif // LR_257
#ifndef RF_257
    #define RF_257
#endif // RF_257
#ifndef DB_257
    #define DB_257
#endif // DB_257
#ifndef ES_258
    #define ES_258
#endif // ES_258
#ifndef FM_258
    #define FM_258
#endif // FM_258
#ifndef MD_258
    #define MD_258
#endif // MD_258
#ifndef NO_258
    #define NO_258
#endif // NO_258
#ifndef LR_258
    #define LR_258
#endif // LR_258
#ifndef RF_258
    #define RF_258
#endif // RF_258
#ifndef DB_258
    #define DB_258
#endif // DB_258
#ifndef SG_258
    #define SG_258
#endif // SG_258
#ifndef FM_259
    #define FM_259
#endif // FM_259
#ifndef MD_259
    #define MD_259
#endif // MD_259
#ifndef NO_259
    #define NO_259
#endif // NO_259
#ifndef LR_259
    #define LR_259
#endif // LR_259
#ifndef RF_259
    #define RF_259
#endif // RF_259
#ifndef DB_259
    #define DB_259
#endif // DB_259
#ifndef SG_259
    #define SG_259
#endif // SG_259
#ifndef FM_260
    #define FM_260
#endif // FM_260
#ifndef MD_260
    #define MD_260
#endif // MD_260
#ifndef NO_260
    #define NO_260
#endif // NO_260
#ifndef LR_260
    #define LR_260
#endif // LR_260
#ifndef RF_260
    #define RF_260
#endif // RF_260
#ifndef DB_260
    #define DB_260
#endif // DB_260
#ifndef SG_260
    #define SG_260
#endif // SG_260
#ifndef BH_260
    #define BH_260
#endif // BH_260
#ifndef MD_261
    #define MD_261
#endif // MD_261
#ifndef NO_261
    #define NO_261
#endif // NO_261
#ifndef LR_261
    #define LR_261
#endif // LR_261
#ifndef RF_261
    #define RF_261
#endif // RF_261
#ifndef DB_261
    #define DB_261
#endif // DB_261
#ifndef SG_261
    #define SG_261
#endif // SG_261
#ifndef BH_261
    #define BH_261
#endif // BH_261
#ifndef MD_262
    #define MD_262
#endif // MD_262
#ifndef NO_262
    #define NO_262
#endif // NO_262
#ifndef LR_262
    #define LR_262
#endif // LR_262
#ifndef RF_262
    #define RF_262
#endif // RF_262
#ifndef DB_262
    #define DB_262
#endif // DB_262
#ifndef SG_262
    #define SG_262
#endif // SG_262
#ifndef BH_262
    #define BH_262
#endif // BH_262
#ifndef NO_263
    #define NO_263
#endif // NO_263
#ifndef LR_263
    #define LR_263
#endif // LR_263
#ifndef RF_263
    #define RF_263
#endif // RF_263
#ifndef DB_263
    #define DB_263
#endif // DB_263
#ifndef SG_263
    #define SG_263
#endif // SG_263
#ifndef BH_263
    #define BH_263
#endif // BH_263
#ifndef HS_263
    #define HS_263
#endif // HS_263
#ifndef NO_264
    #define NO_264
#endif // NO_264
#ifndef LR_264
    #define LR_264
#endif // LR_264
#ifndef RF_264
    #define RF_264
#endif // RF_264
#ifndef DB_264
    #define DB_264
#endif // DB_264
#ifndef SG_264
    #define SG_264
#endif // SG_264
#ifndef BH_264
    #define BH_264
#endif // BH_264
#ifndef HS_264
    #define HS_264
#endif // HS_264
#ifndef LR_265
    #define LR_265
#endif // LR_265
#ifndef RF_265
    #define RF_265
#endif // RF_265
#ifndef DB_265
    #define DB_265
#endif // DB_265
#ifndef SG_265
    #define SG_265
#endif // SG_265
#ifndef BH_265
    #define BH_265
#endif // BH_265
#ifndef HS_265
    #define HS_265
#endif // HS_265
#ifndef MT_265
    #define MT_265
#endif // MT_265
#ifndef LR_266
    #define LR_266
#endif // LR_266
#ifndef RF_266
    #define RF_266
#endif // RF_266
#ifndef DB_266
    #define DB_266
#endif // DB_266
#ifndef SG_266
    #define SG_266
#endif // SG_266
#ifndef BH_266
    #define BH_266
#endif // BH_266
#ifndef HS_266
    #define HS_266
#endif // HS_266
#ifndef MT_266
    #define MT_266
#endif // MT_266
#ifndef RF_267
    #define RF_267
#endif // RF_267
#ifndef DB_267
    #define DB_267
#endif // DB_267
#ifndef SG_267
    #define SG_267
#endif // SG_267
#ifndef BH_267
    #define BH_267
#endif // BH_267
#ifndef HS_267
    #define HS_267
#endif // HS_267
#ifndef MT_267
    #define MT_267
#endif // MT_267
#ifndef DS_267
    #define DS_267
#endif // DS_267
#ifndef RF_268
    #define RF_268
#endif // RF_268
#ifndef DB_268
    #define DB_268
#endif // DB_268
#ifndef SG_268
    #define SG_268
#endif // SG_268
#ifndef BH_268
    #define BH_268
#endif // BH_268
#ifndef HS_268
    #define HS_268
#endif // HS_268
#ifndef MT_268
    #define MT_268
#endif // MT_268
#ifndef DS_268
    #define DS_268
#endif // DS_268
#ifndef DB_269
    #define DB_269
#endif // DB_269
#ifndef SG_269
    #define SG_269
#endif // SG_269
#ifndef BH_269
    #define BH_269
#endif // BH_269
#ifndef HS_269
    #define HS_269
#endif // HS_269
#ifndef MT_269
    #define MT_269
#endif // MT_269
#ifndef DS_269
    #define DS_269
#endif // DS_269
#ifndef DB_270
    #define DB_270
#endif // DB_270
#ifndef SG_270
    #define SG_270
#endif // SG_270
#ifndef BH_270
    #define BH_270
#endif // BH_270
#ifndef HS_270
    #define HS_270
#endif // HS_270
#ifndef MT_270
    #define MT_270
#endif // MT_270
#ifndef DS_270
    #define DS_270
#endif // DS_270
#ifndef SG_271
    #define SG_271
#endif // SG_271
#ifndef BH_271
    #define BH_271
#endif // BH_271
#ifndef HS_271
    #define HS_271
#endif // HS_271
#ifndef MT_271
    #define MT_271
#endif // MT_271
#ifndef DS_271
    #define DS_271
#endif // DS_271
#ifndef SG_272
    #define SG_272
#endif // SG_272
#ifndef BH_272
    #define BH_272
#endif // BH_272
#ifndef HS_272
    #define HS_272
#endif // HS_272
#ifndef MT_272
    #define MT_272
#endif // MT_272
#ifndef DS_272
    #define DS_272
#endif // DS_272
#ifndef RG_272
    #define RG_272
#endif // RG_272
#ifndef SG_273
    #define SG_273
#endif // SG_273
#ifndef BH_273
    #define BH_273
#endif // BH_273
#ifndef HS_273
    #define HS_273
#endif // HS_273
#ifndef MT_273
    #define MT_273
#endif // MT_273
#ifndef DS_273
    #define DS_273
#endif // DS_273
#ifndef RG_273
    #define RG_273
#endif // RG_273
#ifndef BH_274
    #define BH_274
#endif // BH_274
#ifndef HS_274
    #define HS_274
#endif // HS_274
#ifndef MT_274
    #define MT_274
#endif // MT_274
#ifndef DS_274
    #define DS_274
#endif // DS_274
#ifndef RG_274
    #define RG_274
#endif // RG_274
#ifndef BH_275
    #define BH_275
#endif // BH_275
#ifndef HS_275
    #define HS_275
#endif // HS_275
#ifndef MT_275
    #define MT_275
#endif // MT_275
#ifndef DS_275
    #define DS_275
#endif // DS_275
#ifndef RG_275
    #define RG_275
#endif // RG_275
#ifndef BH_276
    #define BH_276
#endif // BH_276
#ifndef HS_276
    #define HS_276
#endif // HS_276
#ifndef MT_276
    #define MT_276
#endif // MT_276
#ifndef DS_276
    #define DS_276
#endif // DS_276
#ifndef RG_276
    #define RG_276
#endif // RG_276
#ifndef CN_276
    #define CN_276
#endif // CN_276
#ifndef BH_277
    #define BH_277
#endif // BH_277
#ifndef HS_277
    #define HS_277
#endif // HS_277
#ifndef MT_277
    #define MT_277
#endif // MT_277
#ifndef DS_277
    #define DS_277
#endif // DS_277
#ifndef RG_277
    #define RG_277
#endif // RG_277
#ifndef CN_277
    #define CN_277
#endif // CN_277
#ifndef BH_278
    #define BH_278
#endif // BH_278
#ifndef HS_278
    #define HS_278
#endif // HS_278
#ifndef MT_278
    #define MT_278
#endif // MT_278
#ifndef DS_278
    #define DS_278
#endif // DS_278
#ifndef RG_278
    #define RG_278
#endif // RG_278
#ifndef CN_278
    #define CN_278
#endif // CN_278
#ifndef NH_278
    #define NH_278
#endif // NH_278
#ifndef HS_279
    #define HS_279
#endif // HS_279
#ifndef MT_279
    #define MT_279
#endif // MT_279
#ifndef DS_279
    #define DS_279
#endif // DS_279
#ifndef RG_279
    #define RG_279
#endif // RG_279
#ifndef CN_279
    #define CN_279
#endif // CN_279
#ifndef NH_279
    #define NH_279
#endif // NH_279
#ifndef HS_280
    #define HS_280
#endif // HS_280
#ifndef MT_280
    #define MT_280
#endif // MT_280
#ifndef DS_280
    #define DS_280
#endif // DS_280
#ifndef RG_280
    #define RG_280
#endif // RG_280
#ifndef CN_280
    #define CN_280
#endif // CN_280
#ifndef NH_280
    #define NH_280
#endif // NH_280
#ifndef MT_281
    #define MT_281
#endif // MT_281
#ifndef DS_281
    #define DS_281
#endif // DS_281
#ifndef RG_281
    #define RG_281
#endif // RG_281
#ifndef CN_281
    #define CN_281
#endif // CN_281
#ifndef NH_281
    #define NH_281
#endif // NH_281
#ifndef MT_282
    #define MT_282
#endif // MT_282
#ifndef DS_282
    #define DS_282
#endif // DS_282
#ifndef RG_282
    #define RG_282
#endif // RG_282
#ifndef CN_282
    #define CN_282
#endif // CN_282
#ifndef NH_282
    #define NH_282
#endif // NH_282
#ifndef DS_283
    #define DS_283
#endif // DS_283
#ifndef RG_283
    #define RG_283
#endif // RG_283
#ifndef CN_283
    #define CN_283
#endif // CN_283
#ifndef NH_283
    #define NH_283
#endif // NH_283
#ifndef DS_284
    #define DS_284
#endif // DS_284
#ifndef RG_284
    #define RG_284
#endif // RG_284
#ifndef CN_284
    #define CN_284
#endif // CN_284
#ifndef NH_284
    #define NH_284
#endif // NH_284
#ifndef FL_284
    #define FL_284
#endif // FL_284
#ifndef RG_285
    #define RG_285
#endif // RG_285
#ifndef CN_285
    #define CN_285
#endif // CN_285
#ifndef NH_285
    #define NH_285
#endif // NH_285
#ifndef FL_285
    #define FL_285
#endif // FL_285
#ifndef RG_286
    #define RG_286
#endif // RG_286
#ifndef CN_286
    #define CN_286
#endif // CN_286
#ifndef NH_286
    #define NH_286
#endif // NH_286
#ifndef FL_286
    #define FL_286
#endif // FL_286
#ifndef CN_287
    #define CN_287
#endif // CN_287
#ifndef NH_287
    #define NH_287
#endif // NH_287
#ifndef FL_287
    #define FL_287
#endif // FL_287
#ifndef MC_287
    #define MC_287
#endif // MC_287
#ifndef CN_288
    #define CN_288
#endif // CN_288
#ifndef NH_288
    #define NH_288
#endif // NH_288
#ifndef FL_288
    #define FL_288
#endif // FL_288
#ifndef MC_288
    #define MC_288
#endif // MC_288
#ifndef NH_289
    #define NH_289
#endif // NH_289
#ifndef FL_289
    #define FL_289
#endif // FL_289
#ifndef MC_289
    #define MC_289
#endif // MC_289
#ifndef LV_289
    #define LV_289
#endif // LV_289
#ifndef NH_290
    #define NH_290
#endif // NH_290
#ifndef FL_290
    #define FL_290
#endif // FL_290
#ifndef MC_290
    #define MC_290
#endif // MC_290
#ifndef LV_290
    #define LV_290
#endif // LV_290
#ifndef FL_291
    #define FL_291
#endif // FL_291
#ifndef MC_291
    #define MC_291
#endif // MC_291
#ifndef LV_291
    #define LV_291
#endif // LV_291
#ifndef TS_291
    #define TS_291
#endif // TS_291
#ifndef MC_292
    #define MC_292
#endif // MC_292
#ifndef LV_292
    #define LV_292
#endif // LV_292
#ifndef TS_292
    #define TS_292
#endif // TS_292
#ifndef LV_293
    #define LV_293
#endif // LV_293
#ifndef TS_293
    #define TS_293
#endif // TS_293
#ifndef OG_293
    #define OG_293
#endif // OG_293
#ifndef TS_294
    #define TS_294
#endif // TS_294
#ifndef OG_294
    #define OG_294
#endif // OG_294
#ifndef OG_295
    #define OG_295
#endif // OG_295
#endif // SPECIES_MASS_DATA_H
