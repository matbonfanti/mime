from math import *


def write_mctdh_opfile(nMonomer, nBath, GSfunc, ESfunc, COUPfunc, RingGS_Morse, RingES_Morse):
    # prepare 1D and 2D functions needed to construct the hamiltonian

    torsion_eq = 45. / 180. * pi

    # 2D function for ground state potential as a function of torsion and stretching
    GS_2D = GSfunc.symbolic_expr().as_coefficients_dict()
    # 1D function for the ground state potential
    GS_1D_x = GSfunc.symbolic_expr().subs(GSfunc.y, torsion_eq).as_coefficients_dict()

    # 2D function for excited state potential as a function of torsion and stretching
    ES_2D = ESfunc.symbolic_expr().as_coefficients_dict()
    # 1D function for the excited state potential as a function of stretching
    ES_1D_x = ESfunc.symbolic_expr().subs(ESfunc.y, torsion_eq).as_coefficients_dict()

    # 2D function for coupling potential as a function of torsion and stretching
    COUP_2D = COUPfunc.symbolic_expr().as_coefficients_dict()
    # 1D function for the coupling potential as a function of stretching
    COUP_1D_x = COUPfunc.symbolic_expr().subs(COUPfunc.y, torsion_eq).as_coefficients_dict()

    # define dictionary for potential printing, GROUND STATE
    GS_1D_x_print = []
    for monom, n in zip(GS_1D_x, range(len(GS_1D_x))):
        GS_1D_x_print.append(["gs_str_a{}".format(n), GS_1D_x[monom], GSfunc.xBasisFmtDict[monom]])
    GS_2D_print = []
    for term, n in zip(GS_2D, range(len(GS_2D))):
        fun_tor, fun_str = term.as_independent(GSfunc.x)
        GS_2D_print.append(
            ["gs_strtor_b{}".format(n), GS_2D[term], GSfunc.xBasisFmtDict[fun_str], GSfunc.yBasisFmtDict[fun_tor]])

    # define dictionary for potential printing, EXCITED STATE
    ES_1D_x_print = []
    for monom, n in zip(ES_1D_x, range(len(ES_1D_x))):
        ES_1D_x_print.append(["es_str_a{}".format(n), ES_1D_x[monom], ESfunc.xBasisFmtDict[monom]])
    ES_2D_print = []
    for term, n in zip(ES_2D, range(len(ES_2D))):
        fun_tor, fun_str = term.as_independent(ESfunc.x)
        ES_2D_print.append(
            ["es_strtor_b{}".format(n), ES_2D[term], ESfunc.xBasisFmtDict[fun_str], ESfunc.yBasisFmtDict[fun_tor]])

    # define dictionary for potential printing, COUPLING
    COUP_1D_x_print = []
    for monom, n in zip(COUP_1D_x, range(len(COUP_1D_x))):
        COUP_1D_x_print.append(["coup_str_a{}".format(n), COUP_1D_x[monom], COUPfunc.xBasisFmtDict[monom]])
    COUP_2D_print = []
    for term, n in zip(COUP_2D, range(len(COUP_2D))):
        fun_tor, fun_str = term.as_independent(COUPfunc.x)
        COUP_2D_print.append(["coup_strtor_b{}".format(n), COUP_2D[term], COUPfunc.xBasisFmtDict[fun_str],
                              COUPfunc.yBasisFmtDict[fun_tor]])

    # define function to print lines for 1D
    def stringfor1Dpotential(factor, nr_dof, potential_print, el=None):
        string = ""
        for t in potential_print:
            if el is None:
                string += "{}*{} |{} {} \n".format(factor, t[0], nr_dof, t[2])
            else:
                string += "{}*{} |1 S{}&{} |{} {} \n".format(factor, t[0], el[0], el[1], nr_dof, t[2])
        return string

    # define function to print lines for 2D
    def stringfor2Dpotential(factor, nr_dof1, nr_dof2, potential_print, el=None):
        string = ""
        for t in potential_print:
            if el is None:
                string += "{}*{} |{} {} |{} {} \n".format(factor, t[0], nr_dof1, t[2], nr_dof2, t[3])
            else:
                string += "{}*{} |1 S{}&{} |{} {} |{} {} \n".format(factor, t[0], el[0], el[1], nr_dof1, t[2],
                                                                    nr_dof2, t[3])
        return string

    f = open(str(nMonomer) + "-mer.op", "w")

    # list of the torsions that are not fixed in the dynamics (nr of the left monomer)
    MovingTorsion = [nMonomer / 2]

    # define names of the variables
    elcoord = "el"
    act_torsions = ["tor" + "%02d%02d" % (i, i + 1) for i in MovingTorsion]
    b_stretch = ["bs" + "%02d%02d" % (i + 1, i + 2) for i in range(nMonomer - 1)]
    internal = ["int" + "%02d" % (i + 1) for i in range(nMonomer)]
    bath_modes = ["bath" + "%02d" % (i + 1) for i in range(nBath * len(MovingTorsion))]

    # define numbers of the variables
    nr_elcoord = 1
    nr_act_torsions = [nr_elcoord + 1 + i for i in range(len(MovingTorsion))]
    nr_b_stretch = [nr_act_torsions[-1] + 1 + i for i in range(len(b_stretch))]
    nr_internal = [nr_b_stretch[-1] + 1 + i for i in range(len(internal))]
    nr_bath_modes = [nr_internal[-1] + 1 + i for i in range(len(bath_modes))]

    f.write("OP_DEFINE-section\ntitle\n1stack of {}mers\nend-title\nend-op_define-section \n""".format(nMonomer))

    f.write("\nPARAMETER-SECTION\n")

    f.write("\n# masses of the torsions\n")
    for label in act_torsions:
        f.write("mass_{} = 1365700.43360989".format(label) + "\n")
    f.write("\n# masses of the bond stretching modes \n")
    for label in b_stretch:
        f.write("mass_{} = 13985.5445411623".format(label) + "\n")
    f.write("\n# masses of the internal modes \n")
    for label in internal:
        f.write("mass_{} = 94846.9644229036".format(label) + "\n")
    f.write("\n# constant for the kinetic coupling between stretching and internal modes \n")
    f.write("gbsint = -0.0000000109\n")

    f.write("\n# coefficients for the ground state potential as a function of the stretching mode \n")
    for term in GS_1D_x_print:
        f.write("{} = {}\n".format(term[0], term[1]))
    f.write("\n# coefficients for the ground state potential as a function of the stretching mode and torsion mode \n")
    for term in GS_2D_print:
        f.write("{} = {}\n".format(term[0], term[1]))

    f.write("\n# coefficients for the excited state potential as a function of the stretching mode \n")
    for term in ES_1D_x_print:
        f.write("{} = {}\n".format(term[0], term[1]))
    f.write("\n# coefficients for the excited state potential as a function of the stretching mode and torsion mode \n")
    for term in ES_2D_print:
        f.write("{} = {}\n".format(term[0], term[1]))

    f.write("\n# coefficients for the coupling potential as a function of the stretching mode \n")
    for term in COUP_1D_x_print:
        f.write("{} = {}\n".format(term[0], term[1]))
    f.write("\n# coefficients for the coupling potential as a function of the stretching mode and torsion mode \n")
    for term in COUP_2D_print:
        f.write("{} = {}\n".format(term[0], term[1]))

    f.write("\n# coefficients for the ground state potential as a function of the ring breathing mode \n")
    f.write("{} = {}\n".format("digs", RingGS_Morse[0]))
    f.write("{} = {}\n".format("alphaigs", RingGS_Morse[1]))
    f.write("{} = {}\n".format("x0igs", RingGS_Morse[2]))
    f.write("{} = {}\n".format("e0igs", RingGS_Morse[3]))
    f.write("\n# coefficients for the excited state potential as a function of the ring breathing mode \n")
    f.write("{} = {}\n".format("dies", RingES_Morse[0]))
    f.write("{} = {}\n".format("alphaies", RingES_Morse[1]))
    f.write("{} = {}\n".format("x0ies", RingES_Morse[2]))
    f.write("{} = {}\n".format("e0ies", RingES_Morse[3]))

    f.write("\n# parameters for the definition of the bath modes \n")
    f.write("mbath = 1.0\n")
    f.write("deltaom = 0.00003\n")
    f.write("cfac = 0.0355222904\n")
    f.write("shi = cfac*cfac/mbath\n")

    f.write("\n# frequencies of the bath modes \n")
    for i in range(len(MovingTorsion)):
        for j in range(nBath):
            f.write("om{:02d} = {}*deltaom\n".format(i * nBath + j + 1, float(j + 1)))
    f.write("\n# masses of the bath modes \n")
    for label in bath_modes:
        f.write("mass_{} = mbath".format(label) + "\n")
    f.write("\n# quadratic potential terms of the bath modes \n")
    for i in range(len(MovingTorsion) * nBath):
        f.write("k{:02d} = 0.5*mbath*om{:02d}*om{:02d}\n".format(i + 1, i + 1, i + 1))
    f.write("\n# coupling coefficients to the system \n")
    for i in range(len(MovingTorsion) * nBath):
        f.write("c{:02d} = cfac*om{:02d}\n".format(i + 1, i + 1))

    f.write("\nend-parameter-section\n")

    f.write("""\nLABELS-SECTION

# morse potentials for the ground and excited state momomer V as a function of the ring breathing internal mode
intgs = morse1[digs,alphaigs,x0igs,e0igs]
intes = morse1[dies,alphaies,x0ies,e0ies]

end-labels-section
""")

    f.write("\nHAMILTONIAN-SECTION\n")
    f.write("\n# electronic mode\n")
    f.write("modes | {}".format(elcoord) + "\n")

    def write_formatted_names_of_modes(listofnames, size):
        for sublist in [listofnames[i:i + size] for i in range(0, len(listofnames), size)]:
            f.write("modes " + " ".join(["| {:^7}".format(s) for s in sublist]) + "\n")

    f.write("# name for the central torsion mode\n")
    write_formatted_names_of_modes(act_torsions, 10)
    f.write("# names for the bond stretching modes\n")
    write_formatted_names_of_modes(b_stretch, 10)
    f.write("# names for the internal monomer modes\n")
    write_formatted_names_of_modes(internal, 10)
    f.write("# names for the bath modes coupled to central torsion\n")
    write_formatted_names_of_modes(bath_modes, 10)

    f.write("\n# kinetic coupling between stretching and internal modes\n")
    for i in range(nMonomer - 1):
        f.write("-gbsint |{} dq |{} dq\n".format(nr_b_stretch[i], nr_internal[i]))
        f.write("-gbsint |{} dq |{} dq\n".format(nr_b_stretch[i], nr_internal[i + 1]))

    f.write("\n# Kinetic operator for the torsion\n")
    for i in nr_act_torsions:
        f.write("1.0 |{} KE\n".format(i))
    f.write("\n# Kinetic operator for the stretching coordinates\n")
    for i in nr_b_stretch:
        f.write("1.0 |{} KE\n".format(i))
    f.write("\n# Kinetic operator for the internal monomer coordinates\n")
    for i in nr_internal:
        f.write("1.0 |{} KE\n".format(i))

    f.write("\n# Ground state potential\n")
    f.write("# potential of the internal modes\n")
    for i in range(nMonomer):
        f.write("1.0 |{} {}\n".format(nr_internal[i], "intgs"))
    f.write("# potential of the inter-monomer stretching and active torsion modes\n")
    for i in range(nMonomer - 1):
        if i + 1 in MovingTorsion:
            tor_i = MovingTorsion.index(i + 1)
            f.write(stringfor2Dpotential(2.0, nr_b_stretch[i], nr_act_torsions[tor_i], GS_2D_print))
        else:
            f.write(stringfor1Dpotential(2.0, nr_b_stretch[i], GS_1D_x_print))

    for i in range(nMonomer):
        f.write("\n# Diagonal excited state potential, exciton on the monomer {}\n".format(i + 1))
        f.write("# potential of the internal modes\n")
        f.write(" 1.0 |1 S{}&{} |{} {} \n".format(i + 1, i + 1, nr_internal[i], "intes"))
        f.write("-1.0 |1 S{}&{} |{} {} \n".format(i + 1, i + 1, nr_internal[i], "intgs"))
        if i != 0:  # this is not the first monomer, left torsion and stretching do exist
            f.write("# potential of the left intermonomer modes\n")
            if i in MovingTorsion:
                tor_i = MovingTorsion.index(i)
                f.write(
                    stringfor2Dpotential(1.0, nr_b_stretch[i - 1], nr_act_torsions[tor_i], ES_2D_print, [i + 1, i + 1]))
                f.write(stringfor2Dpotential(-1.0, nr_b_stretch[i - 1], nr_act_torsions[tor_i], GS_2D_print,
                                             [i + 1, i + 1]))
            else:
                f.write(stringfor1Dpotential(1.0, nr_b_stretch[i - 1], ES_1D_x_print, [i + 1, i + 1]))
                f.write(stringfor1Dpotential(-1.0, nr_b_stretch[i - 1], GS_1D_x_print, [i + 1, i + 1]))
        if i != nMonomer - 1:  # this is not the last monomer, right torsion and stretching do exist
            f.write("# potential of the right intermonomer modes\n")
            if i + 1 in MovingTorsion:
                tor_i = MovingTorsion.index(i + 1)
                f.write(stringfor2Dpotential(1.0, nr_b_stretch[i], nr_act_torsions[tor_i], ES_2D_print, [i + 1, i + 1]))
                f.write(
                    stringfor2Dpotential(-1.0, nr_b_stretch[i], nr_act_torsions[tor_i], GS_2D_print, [i + 1, i + 1]))
            else:
                f.write(stringfor1Dpotential(1.0, nr_b_stretch[i], ES_1D_x_print, [i + 1, i + 1]))
                f.write(stringfor1Dpotential(-1.0, nr_b_stretch[i], GS_1D_x_print, [i + 1, i + 1]))

    for i in range(nMonomer - 1):
        f.write("\n# Coupling potential between excitons on monomer {} and {}\n".format(i + 1, i + 2))
        f.write("# potential of the intermonomer modes\n")
        if i + 1 in MovingTorsion:
            tor_i = MovingTorsion.index(i + 1)
            f.write(stringfor2Dpotential(1.0, nr_b_stretch[i], nr_act_torsions[tor_i], COUP_2D_print, [i + 1, i + 2]))
        else:
            f.write(stringfor1Dpotential(1.0, nr_b_stretch[i], COUP_1D_x_print, [i + 1, i + 2]))

    f.write("\n# Kinetic energy of the bath coupled to the torsion of mode {}\n".format(1))
    for i in nr_bath_modes:
        f.write("1.0 |{} KE\n".format(i))

    f.write("\n# Potential of the bath coupled to the torsion of mode {}\n".format(1))
    for i, nr in zip(range(len(MovingTorsion) * nBath), nr_bath_modes):
        f.write("k{:02d} |{} q^2\n".format(i + 1, nr))
    for i in range(len(MovingTorsion)):
        for j in range(nBath):
            f.write(
                "-c{:02d} |{} q |{} q \n".format(i * nBath + j + 1, nr_act_torsions[i], nr_bath_modes[i * nBath + j]))
    for i in range(len(MovingTorsion)):
        f.write("20.0*0.5*shi   |{} q^2 \n".format(nr_act_torsions[i]))

    f.write("\nend-hamiltonian-section\n")
    f.write("\nend-operator\n")

    f.close()
