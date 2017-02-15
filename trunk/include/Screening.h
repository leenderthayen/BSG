namespace screening {

//================================================================================
// Salvat's parameters for the screened potential

inline void PotParam(int Zloc, std::vector<double> &Aby, std::vector<double> &Bby) {
  std::vector<double> Aby_loc, Bby_loc;

  Aby.clear();
  Bby.clear();

  double b = 0.;
  if (Zloc < 0) Zloc = -Zloc;  // for Z<0 if beta + transition

  if (Zloc > 92)  // Moliere's potential
  {
    b = 0.88534 * std::pow(Zloc * 1., -1. / 3.);
    Aby_loc.push_back(0.1);
    Aby_loc.push_back(0.55);
    Aby_loc.push_back(0.35);
    Bby_loc.push_back(6.0 / b);
    Bby_loc.push_back(1.2 / b);
    Bby_loc.push_back(0.3 / b);
  }

  else {
    // Screening for H1
    if (Zloc == 1) {
      Aby_loc.push_back(-184.39);
      Aby_loc.push_back(185.39);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(2.0027);
      Bby_loc.push_back(1.9973);
      Bby_loc.push_back(0.0000);
    }

    // Screening for He2
    if (Zloc == 2) {
      Aby_loc.push_back(-0.2259);
      Aby_loc.push_back(1.2259);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(5.5272);
      Bby_loc.push_back(2.3992);
      Bby_loc.push_back(0.0000);
    }

    // Screening for Li3
    if (Zloc == 3) {
      Aby_loc.push_back(0.6045);
      Aby_loc.push_back(0.3955);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(2.8174);
      Bby_loc.push_back(0.6625);
      Bby_loc.push_back(0.0000);
    }

    // Screening for Be4
    if (Zloc == 4) {
      Aby_loc.push_back(0.3278);
      Aby_loc.push_back(0.6722);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(04.5430);
      Bby_loc.push_back(0.9852);
      Bby_loc.push_back(0.0000);
    }

    // Screening for B5
    if (Zloc == 5) {
      Aby_loc.push_back(0.2327);
      Aby_loc.push_back(0.7673);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(5.9900);
      Bby_loc.push_back(1.2135);
      Bby_loc.push_back(0.0000);
    }

    // Screening for C6
    if (Zloc == 6) {
      Aby_loc.push_back(0.1537);
      Aby_loc.push_back(0.8463);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(8.0404);
      Bby_loc.push_back(1.4913);
      Bby_loc.push_back(0.0000);
    }

    // Screening for N7
    if (Zloc == 7) {
      Aby_loc.push_back(0.0996);
      Aby_loc.push_back(0.9004);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(10.812);
      Bby_loc.push_back(1.7687);
      Bby_loc.push_back(0.0000);
    }

    // Screening for O8
    if (Zloc == 8) {
      Aby_loc.push_back(0.0625);
      Aby_loc.push_back(0.9375);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(14.823);
      Bby_loc.push_back(2.0403);
      Bby_loc.push_back(0.0000);
    }

    // Screening for F9
    if (Zloc == 9) {
      Aby_loc.push_back(0.0368);
      Aby_loc.push_back(0.9632);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(21.400);
      Bby_loc.push_back(2.3060);
      Bby_loc.push_back(0.0000);
    }

    // Screening for Ne10
    if (Zloc == 10) {
      Aby_loc.push_back(0.0188);
      Aby_loc.push_back(0.9812);
      Aby_loc.push_back(1 - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(34.999);
      Bby_loc.push_back(2.5662);
      Bby_loc.push_back(0.0000);
    }

    // Screening for Na11
    if (Zloc == 11) {
      Aby_loc.push_back(0.7444);
      Aby_loc.push_back(0.2556);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(4.1205);
      Bby_loc.push_back(0.8718);
      Bby_loc.push_back(0.0000);
    }

    // Screening for Mg12
    if (Zloc == 12) {
      Aby_loc.push_back(0.6423);
      Aby_loc.push_back(0.3577);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(4.7266);
      Bby_loc.push_back(1.0025);
      Bby_loc.push_back(0.0000);
    }

    // Screening for Al13
    if (Zloc == 13) {
      Aby_loc.push_back(0.6002);
      Aby_loc.push_back(0.3998);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(5.1405);
      Bby_loc.push_back(1.0153);
      Bby_loc.push_back(0.0000);
    }

    // Screening for Si14
    if (Zloc == 14) {
      Aby_loc.push_back(0.5160);
      Aby_loc.push_back(0.4840);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(5.8492);
      Bby_loc.push_back(1.1732);
      Bby_loc.push_back(0.0000);
    }

    // Screening for P15
    if (Zloc == 15) {
      Aby_loc.push_back(0.4387);
      Aby_loc.push_back(0.5613);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(6.6707);
      Bby_loc.push_back(1.3410);
      Bby_loc.push_back(0.0000);
    }

    // Screening for S16
    if (Zloc == 16) {
      Aby_loc.push_back(0.5459);
      Aby_loc.push_back(-0.5333);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(6.3703);
      Bby_loc.push_back(2.5517);
      Bby_loc.push_back(1.6753);
    }

    // Screening for Cl17
    if (Zloc == 17) {
      Aby_loc.push_back(0.7249);
      Aby_loc.push_back(-0.7548);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(6.2118);
      Bby_loc.push_back(3.3883);
      Bby_loc.push_back(1.8596);
    }

    // Screening for Ar18
    if (Zloc == 18) {
      Aby_loc.push_back(2.1912);
      Aby_loc.push_back(-2.2852);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(5.5470);
      Bby_loc.push_back(4.5687);
      Bby_loc.push_back(2.0446);
    }

    // Screening for K19
    if (Zloc == 19) {
      Aby_loc.push_back(0.0486);
      Aby_loc.push_back(0.7759);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(30.260);
      Bby_loc.push_back(3.1243);
      Bby_loc.push_back(0.7326);
    }

    // Screening for Ca20
    if (Zloc == 20) {
      Aby_loc.push_back(0.5800);
      Aby_loc.push_back(0.4200);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(6.3218);
      Bby_loc.push_back(1.0094);
      Bby_loc.push_back(0.0000);
    }

    // Screening for Sc21
    if (Zloc == 21) {
      Aby_loc.push_back(0.5543);
      Aby_loc.push_back(0.4457);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(6.6328);
      Bby_loc.push_back(1.1023);
      Bby_loc.push_back(0.0000);
    }

    // Screening for Ti22
    if (Zloc == 22) {
      Aby_loc.push_back(0.0112);
      Aby_loc.push_back(0.6832);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(99.757);
      Bby_loc.push_back(4.1286);
      Bby_loc.push_back(1.0090);
    }

    // Screening for V23
    if (Zloc == 23) {
      Aby_loc.push_back(0.0318);
      Aby_loc.push_back(0.6753);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(42.533);
      Bby_loc.push_back(3.9404);
      Bby_loc.push_back(1.0533);
    }

    // Screening for Cr24
    if (Zloc == 24) {
      Aby_loc.push_back(0.1075);
      Aby_loc.push_back(0.7162);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(18.959);
      Bby_loc.push_back(3.0638);
      Bby_loc.push_back(1.0014);
    }

    // Screening for Mn25
    if (Zloc == 25) {
      Aby_loc.push_back(0.0498);
      Aby_loc.push_back(0.6866);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(31.864);
      Bby_loc.push_back(3.7811);
      Bby_loc.push_back(1.1279);
    }

    // Screening for Fe26
    if (Zloc == 26) {
      Aby_loc.push_back(0.0512);
      Aby_loc.push_back(0.6995);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(31.825);
      Bby_loc.push_back(3.7716);
      Bby_loc.push_back(1.1606);
    }

    // Screening for Co27
    if (Zloc == 27) {
      Aby_loc.push_back(0.0500);
      Aby_loc.push_back(0.7142);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(32.915);
      Bby_loc.push_back(3.7908);
      Bby_loc.push_back(1.1915);
    }
    // Screening for Ni28
    if (Zloc == 28) {
      Aby_loc.push_back(0.0474);
      Aby_loc.push_back(0.7294);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(34.758);
      Bby_loc.push_back(3.8299);
      Bby_loc.push_back(1.2209);
    }

    // Screening for Cu29
    if (Zloc == 29) {
      Aby_loc.push_back(0.0771);
      Aby_loc.push_back(0.7951);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(25.326);
      Bby_loc.push_back(3.3928);
      Bby_loc.push_back(1.1426);
    }

    // Screening for Zlocn30
    if (Zloc == 30) {
      Aby_loc.push_back(0.0400);
      Aby_loc.push_back(0.7590);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(40.343);
      Bby_loc.push_back(3.9465);
      Bby_loc.push_back(1.2759);
    }
    // Screening for Ga31
    if (Zloc == 31) {
      Aby_loc.push_back(0.1083);
      Aby_loc.push_back(0.7489);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(20.192);
      Bby_loc.push_back(3.4733);
      Bby_loc.push_back(1.0064);
    }

    // Screening for Ge32
    if (Zloc == 32) {
      Aby_loc.push_back(0.0610);
      Aby_loc.push_back(0.7157);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(29.200);
      Bby_loc.push_back(4.1252);
      Bby_loc.push_back(1.1845);
    }

    // Screening for As33
    if (Zloc == 33) {
      Aby_loc.push_back(0.0212);
      Aby_loc.push_back(0.6709);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(62.487);
      Bby_loc.push_back(4.9502);
      Bby_loc.push_back(1.3582);
    }

    // Screening for Se34
    if (Zloc == 34) {
      Aby_loc.push_back(0.4836);
      Aby_loc.push_back(0.5164);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(8.7824);
      Bby_loc.push_back(1.6967);
      Bby_loc.push_back(0.0000);
    }

    // Screening for Br35
    if (Zloc == 35) {
      Aby_loc.push_back(0.4504);
      Aby_loc.push_back(0.5496);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(9.3348);
      Bby_loc.push_back(1.7900);
      Bby_loc.push_back(0.0000);
    }

    // Screening for Kr36
    if (Zloc == 36) {
      Aby_loc.push_back(0.4190);
      Aby_loc.push_back(0.5810);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(09.9142);
      Bby_loc.push_back(1.8835);
      Bby_loc.push_back(0.0000);
    }

    // Screening for Rb37
    if (Zloc == 37) {
      Aby_loc.push_back(0.1734);
      Aby_loc.push_back(0.7253);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(17.166);
      Bby_loc.push_back(3.1103);
      Bby_loc.push_back(0.7177);
    }

    // Screening for Sr38
    if (Zloc == 38) {
      Aby_loc.push_back(0.0336);
      Aby_loc.push_back(0.7816);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(55.208);
      Bby_loc.push_back(4.2842);
      Bby_loc.push_back(0.8578);
    }

    // Screening for Y39
    if (Zloc == 39) {
      Aby_loc.push_back(0.0689);
      Aby_loc.push_back(0.7202);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(31.366);
      Bby_loc.push_back(4.2412);
      Bby_loc.push_back(0.9472);
    }

    // Screening for Zlocr40
    if (Zloc == 40) {
      Aby_loc.push_back(0.1176);
      Aby_loc.push_back(0.6581);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(22.054);
      Bby_loc.push_back(4.0325);
      Bby_loc.push_back(1.0181);
    }

    // Screening for Nb41
    if (Zloc == 41) {
      Aby_loc.push_back(0.2257);
      Aby_loc.push_back(0.5821);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(14.240);
      Bby_loc.push_back(2.9702);
      Bby_loc.push_back(1.0170);
    }

    // Screening for Mo42
    if (Zloc == 42) {
      Aby_loc.push_back(0.2693);
      Aby_loc.push_back(0.5763);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(14.044);
      Bby_loc.push_back(2.8611);
      Bby_loc.push_back(1.0591);
    }

    // Screening for Tc43
    if (Zloc == 43) {
      Aby_loc.push_back(0.2201);
      Aby_loc.push_back(0.5618);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(15.918);
      Bby_loc.push_back(3.3672);
      Bby_loc.push_back(1.1548);
    }

    // Screening for Ru44
    if (Zloc == 44) {
      Aby_loc.push_back(0.2751);
      Aby_loc.push_back(0.5943);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(14.314);
      Bby_loc.push_back(2.7370);
      Bby_loc.push_back(1.1092);
    }

    // Screening for Rh45
    if (Zloc == 45) {
      Aby_loc.push_back(0.2711);
      Aby_loc.push_back(0.6119);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(14.654);
      Bby_loc.push_back(2.7183);
      Bby_loc.push_back(1.1234);
    }

    // Screening for Pd46
    if (Zloc == 46) {
      Aby_loc.push_back(0.2784);
      Aby_loc.push_back(0.6067);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(14.645);
      Bby_loc.push_back(2.6155);
      Bby_loc.push_back(1.4318);
    }

    // Screening for Ag47
    if (Zloc == 47) {
      Aby_loc.push_back(0.2562);
      Aby_loc.push_back(0.6505);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(15.5880);
      Bby_loc.push_back(2.7412);
      Bby_loc.push_back(1.1408);
    }
    // Screening for Cd48

    if (Zloc == 48) {
      Aby_loc.push_back(0.2271);
      Aby_loc.push_back(0.6155);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(16.914);
      Bby_loc.push_back(3.0841);
      Bby_loc.push_back(1.2619);
    }

    // Screening for In49
    if (Zloc == 49) {
      Aby_loc.push_back(0.2492);
      Aby_loc.push_back(0.6440);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(16.155);
      Bby_loc.push_back(2.8819);
      Bby_loc.push_back(0.9942);
    }

    // Screening for Sn50
    if (Zloc == 50) {
      Aby_loc.push_back(0.2153);
      Aby_loc.push_back(0.6115);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(17.7930);
      Bby_loc.push_back(3.2937);
      Bby_loc.push_back(1.1478);
    }

    // Screening for Sb51
    if (Zloc == 51) {
      Aby_loc.push_back(0.1806);
      Aby_loc.push_back(0.5767);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(19.875);
      Bby_loc.push_back(3.8092);
      Bby_loc.push_back(1.2829);
    }

    // Screening for Te52
    if (Zloc == 52) {
      Aby_loc.push_back(0.1308);
      Aby_loc.push_back(0.5504);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(24.154);
      Bby_loc.push_back(4.6119);
      Bby_loc.push_back(1.4195);
    }

    // Screening for I53
    if (Zloc == 53) {
      Aby_loc.push_back(0.0588);
      Aby_loc.push_back(0.5482);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(39.996);
      Bby_loc.push_back(5.9132);
      Bby_loc.push_back(1.5471);
    }

    // Screening for Xe54
    if (Zloc == 54) {
      Aby_loc.push_back(0.4451);
      Aby_loc.push_back(0.5549);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(11.8050);
      Bby_loc.push_back(1.7967);
      Bby_loc.push_back(0.0000);
    }

    // Screening for Cs55
    if (Zloc == 55) {
      Aby_loc.push_back(0.2708);
      Aby_loc.push_back(0.6524);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(16.591);
      Bby_loc.push_back(2.6964);
      Bby_loc.push_back(0.6814);
    }

    // Screening for Ba56
    if (Zloc == 56) {
      Aby_loc.push_back(0.1728);
      Aby_loc.push_back(0.6845);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(22.397);
      Bby_loc.push_back(3.4595);
      Bby_loc.push_back(0.8073);
    }

    // Screening for La57
    if (Zloc == 57) {
      Aby_loc.push_back(0.1947);
      Aby_loc.push_back(0.6384);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(20.764);
      Bby_loc.push_back(3.4657);
      Bby_loc.push_back(0.8911);
    }

    // Screening for Ce58
    if (Zloc == 58) {
      Aby_loc.push_back(0.1913);
      Aby_loc.push_back(0.6467);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(21.235);
      Bby_loc.push_back(3.4819);
      Bby_loc.push_back(0.9011);
    }

    // Screening for Pr59
    if (Zloc == 59) {
      Aby_loc.push_back(0.1868);
      Aby_loc.push_back(0.6558);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(21.803);
      Bby_loc.push_back(3.5098);
      Bby_loc.push_back(0.9106);
    }

    // Screening for Nd60
    if (Zloc == 60) {
      Aby_loc.push_back(0.1665);
      Aby_loc.push_back(0.7057);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(23.949);
      Bby_loc.push_back(3.5199);
      Bby_loc.push_back(0.8486);
    }

    // Screening for Pm61
    if (Zloc == 61) {
      Aby_loc.push_back(0.1624);
      Aby_loc.push_back(0.7133);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(24.598);
      Bby_loc.push_back(3.5560);
      Bby_loc.push_back(0.8569);
    }

    // Screening for Sm62
    if (Zloc == 62) {
      Aby_loc.push_back(0.1580);
      Aby_loc.push_back(0.7210);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(25.297);
      Bby_loc.push_back(3.5963);
      Bby_loc.push_back(0.8650);
    }

    // Screening for Eu63
    if (Zloc == 63) {
      Aby_loc.push_back(0.1538);
      Aby_loc.push_back(0.7284);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(26.017);
      Bby_loc.push_back(3.6383);
      Bby_loc.push_back(0.8731);
    }

    // Screening for Gd64
    if (Zloc == 64) {
      Aby_loc.push_back(0.1587);
      Aby_loc.push_back(0.7024);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(25.497);
      Bby_loc.push_back(3.7364);
      Bby_loc.push_back(0.9550);
    }

    // Screening for Tb65
    if (Zloc == 65) {
      Aby_loc.push_back(0.1453);
      Aby_loc.push_back(0.7426);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(27.547);
      Bby_loc.push_back(3.7288);
      Bby_loc.push_back(0.8890);
    }

    // Screening for Dy66
    if (Zloc == 66) {
      Aby_loc.push_back(0.1413);
      Aby_loc.push_back(0.7494);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(28.346);
      Bby_loc.push_back(3.7763);
      Bby_loc.push_back(0.8969);
    }

    // Screening for Ho67
    if (Zloc == 67) {
      Aby_loc.push_back(0.1374);
      Aby_loc.push_back(0.7558);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(29.160);
      Bby_loc.push_back(3.8244);
      Bby_loc.push_back(0.9048);
    }

    // Screening for Er68
    if (Zloc == 68) {
      Aby_loc.push_back(0.1336);
      Aby_loc.push_back(0.7619);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(29.990);
      Bby_loc.push_back(3.8734);
      Bby_loc.push_back(0.9128);
    }

    // Screening for Tm69
    if (Zloc == 69) {
      Aby_loc.push_back(0.1299);
      Aby_loc.push_back(0.7680);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(30.835);
      Bby_loc.push_back(3.9233);
      Bby_loc.push_back(0.9203);
    }

    // Screening for Yb70
    if (Zloc == 70) {
      Aby_loc.push_back(0.1267);
      Aby_loc.push_back(0.7734);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(31.681);
      Bby_loc.push_back(3.9727);
      Bby_loc.push_back(0.9288);
    }

    // Screening for Lu71
    if (Zloc == 71) {
      Aby_loc.push_back(0.1288);
      Aby_loc.push_back(0.7528);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(31.353);
      Bby_loc.push_back(4.0904);
      Bby_loc.push_back(1.0072);
    }

    // Screening for Hf72
    if (Zloc == 72) {
      Aby_loc.push_back(0.1303);
      Aby_loc.push_back(0.7324);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(31.217);
      Bby_loc.push_back(4.2049);
      Bby_loc.push_back(1.0946);
    }

    // Screening for Ta73
    if (Zloc == 73) {
      Aby_loc.push_back(0.1384);
      Aby_loc.push_back(0.7096);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(30.077);
      Bby_loc.push_back(4.2492);
      Bby_loc.push_back(1.1697);
    }

    // Screening for W74
    if (Zloc == 74) {
      Aby_loc.push_back(0.1500);
      Aby_loc.push_back(0.6871);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(28.630);
      Bby_loc.push_back(4.2426);
      Bby_loc.push_back(1.2340);
    }

    // Screening for Re75
    if (Zloc == 75) {
      Aby_loc.push_back(0.1608);
      Aby_loc.push_back(0.6659);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(27.568);
      Bby_loc.push_back(4.2341);
      Bby_loc.push_back(1.2970);
    }

    // Screening for Os76
    if (Zloc == 76) {
      Aby_loc.push_back(0.1722);
      Aby_loc.push_back(0.6468);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(26.586);
      Bby_loc.push_back(4.1999);
      Bby_loc.push_back(1.3535);
    }

    // Screening for Ir77
    if (Zloc == 77) {
      Aby_loc.push_back(0.1834);
      Aby_loc.push_back(0.6306);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(25.734);
      Bby_loc.push_back(4.1462);
      Bby_loc.push_back(1.4037);
    }

    // Screening for Pt78
    if (Zloc == 78) {
      Aby_loc.push_back(0.2230);
      Aby_loc.push_back(0.6176);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(22.994);
      Bby_loc.push_back(3.7346);
      Bby_loc.push_back(1.4428);
    }

    // Screening for Au79
    if (Zloc == 79) {
      Aby_loc.push_back(0.2289);
      Aby_loc.push_back(0.6114);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(22.864);
      Bby_loc.push_back(3.6914);
      Bby_loc.push_back(1.4886);
    }

    // Screening for Hg80
    if (Zloc == 80) {
      Aby_loc.push_back(0.2098);
      Aby_loc.push_back(0.6004);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(24.408);
      Bby_loc.push_back(3.9643);
      Bby_loc.push_back(1.5343);
    }

    // Screening for Tl81
    if (Zloc == 81) {
      Aby_loc.push_back(0.2708);
      Aby_loc.push_back(0.6428);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(20.941);
      Bby_loc.push_back(3.2456);
      Bby_loc.push_back(1.1121);
    }

    // Screening for Pb82
    if (Zloc == 82) {
      Aby_loc.push_back(0.2380);
      Aby_loc.push_back(0.6308);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(22.987);
      Bby_loc.push_back(3.6217);
      Bby_loc.push_back(1.2373);
    }

    // Screening for Bi83
    if (Zloc == 83) {
      Aby_loc.push_back(0.2288);
      Aby_loc.push_back(0.6220);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(23.792);
      Bby_loc.push_back(3.7796);
      Bby_loc.push_back(1.2534);
    }

    // Screening for Po84
    if (Zloc == 84) {
      Aby_loc.push_back(0.1941);
      Aby_loc.push_back(0.6105);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(26.695);
      Bby_loc.push_back(4.2582);
      Bby_loc.push_back(1.3577);
    }

    // Screening for At85
    if (Zloc == 85) {
      Aby_loc.push_back(0.1500);
      Aby_loc.push_back(0.6031);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(31.840);
      Bby_loc.push_back(4.9285);
      Bby_loc.push_back(1.4683);
    }

    // Screening for Rn86
    if (Zloc == 86) {
      Aby_loc.push_back(0.0955);
      Aby_loc.push_back(0.6060);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(43.489);
      Bby_loc.push_back(5.8520);
      Bby_loc.push_back(1.5736);
    }

    // Screening for Fr87
    if (Zloc == 87) {
      Aby_loc.push_back(0.3192);
      Aby_loc.push_back(0.6233);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(20.015);
      Bby_loc.push_back(2.9091);
      Bby_loc.push_back(0.7207);
    }

    // Screening for Ra88
    if (Zloc == 88) {
      Aby_loc.push_back(0.2404);
      Aby_loc.push_back(0.6567);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(24.501);
      Bby_loc.push_back(3.5524);
      Bby_loc.push_back(0.8376);
    }

    // Screening for Ac89
    if (Zloc == 89) {
      Aby_loc.push_back(0.2266);
      Aby_loc.push_back(0.6422);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(25.684);
      Bby_loc.push_back(3.7922);
      Bby_loc.push_back(0.9335);
    }

    // Screening for Th90
    if (Zloc == 90) {
      Aby_loc.push_back(0.2176);
      Aby_loc.push_back(0.6240);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(26.554);
      Bby_loc.push_back(4.0044);
      Bby_loc.push_back(1.0238);
    }

    // Screening for Pa91
    if (Zloc == 91) {
      Aby_loc.push_back(0.2413);
      Aby_loc.push_back(0.6304);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(25.193);
      Bby_loc.push_back(3.6780);
      Bby_loc.push_back(0.9699);
    }

    // Screening for U92
    if (Zloc == 92) {
      Aby_loc.push_back(0.2448);
      Aby_loc.push_back(0.6298);
      Aby_loc.push_back(1. - (Aby_loc[0] + Aby_loc[1]));
      Bby_loc.push_back(25.252);
      Bby_loc.push_back(3.6397);
      Bby_loc.push_back(0.9825);
    }
  }

  // Parameters are in atomic units, thus conversion in natural units
  for (int i = 0; i < (int)Bby_loc.size(); i++) Bby_loc[i] = Bby_loc[i] * alpha;

  // Parameters are put in the global std::vectors for the required nucleus
  for (int i = 0; i < (int)Aby_loc.size(); i++) {
    Aby.push_back(Aby_loc[i]);
    Bby.push_back(Bby_loc[i]);
  }

  Aby_loc.clear();
  Bby_loc.clear();
}
}
