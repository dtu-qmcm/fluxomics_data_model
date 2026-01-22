clear functions

r = reaction({...
...
... % Glycolysis and OPP
... % ---------------------------------------------------------------------
      'G6P (abcdef) <-> F6P (abcdef)';...                 % pgi    R1
      'G6P (abcdef) -> Ru5P (bcdef) + CO2 (a)';...        % g6pdh  R2
      'F6P (abcdef) <-> FBP (abcdef)';...                 % pfk    R3
      'FBP (abcdef) <-> DHAP (cba) + GAP (def)';...       % fba    R4
      'DHAP (abc) <-> GAP (abc)';...                      % tpi    R5
      'GAP (abc) <-> G3P (abc)';...                       % gapdh  R6
      'G3P (abc) <-> G2P (abc)';...                       % gpm    R7
      'G2P (abc) <-> PEP (abc)';...                       % eno    R8
      'PEP (abc) <-> Pyr (abc)';...                       % pk     R9
...
... % CBB cycle
... % ---------------------------------------------------------------------
      'Ru5P (abcde) <-> X5P (abcde)';...                  % rpe    R10
      'Ru5P (abcde) <-> R5P (abcde)';...                  % rpi    R11
      'Ru5P (abcde) -> RuBP (abcde)';...                   % prk    R12
      'Ru5P (abcde) + CO2(f) -> G3P (cde) + G3P (fba)';...% rbc1   R13
      'X5P (abcde) <-> TK (ab) + GAP (cde)';...           % tkt1   R14
      'F6P (abcdef) <-> TK (ab) + E4P (cdef)';...         % tkt2   R15
      'S7P (abcdefg) <-> TK (ab) + R5P (cdefg)';...       % tkt3   R16
      'F6P (abcdef) <-> TA (abc) + GAP (def)';...         % tal1   R17
      'S7P (abcdefg) <-> TA (abc) + E4P (defg)';...       % tal2   R18
      'DHAP (abc) + E4P (defg) -> SBP (cbadefg)';...      % sba    R19
      'SBP (abcdefg) -> S7P (abcdefg)';...                % sbp    R20
...
... % TCA cycle
... % ---------------------------------------------------------------------
      'Pyr (abc) -> AcCoA (bc) + CO2 (a)';...             % pdh    R21
      'OAA (abcd) + AcCoA (ef) -> Cit (dcbfea)';...       % cs     R22
      'Cit (abcdef) <-> ICit (abcdef)';...                % can    R23
      'ICit (abcdef) <-> AKG (abcde) + CO2 (f)';...       % icd    R24
      'Suc (abcd) <-> Fum (abcd)';...                     % sdh    R25
      'Fum (abcd) <-> Mal (abcd)';...                     % fum    R26
      'Mal (abcd) <-> OAA (abcd)';...                     % mdh    R27
...
... % Glyoxylate Shunt
... % ---------------------------------------------------------------------
      'ICit (abcdef) -> Glx (ed) + Suc (abcf)';...        % icl    R28
      'Glx (ab) + AcCoA (cd) <-> Mal (abdc)';...          % ms     R29
...
... % Amphibolic reactions
... % ---------------------------------------------------------------------
      'Mal (abcd) -> Pyr (abc) + CO2 (d)';...             % me     R30
      'PEP (abc) + CO2 (d) -> OAA (abcd)';...             % ppc    R31
...
... % Photorespiration
... % ---------------------------------------------------------------------
      'RuBP (abcde) -> G3P (cde) + PG (ba)';...           % rbc2   R32
      'PG (ab) -> Gc (ab)';...                            % pgp    R33
      'Gc (ab) -> Glx (ab)';...                           % gld    R34
      'Glx (ab) + Glx (cd) -> Ga (abc) + CO2 (d)';...     % gt     R35
      'Ga (abc) <-> G3P (abc)';...                        % glyk   R36
...
... % Transport
... % ---------------------------------------------------------------------
      'CO2.ex (a) -> CO2 (a)';...                        % co2in   R37
...
... % Biomass formation
... % ---------------------------------------------------------------------
    %  '0.715 R5P + 3.624 AcCoA + 1.191 G6P + 0.501 E4P + 1.205 G3P + 1.002 PEP + 1.197 Pyr + 2.039 OAA+  1.233 AKG + 0.133 GAP -> Biomass + 0.683 Fum + 1.017 CO2';... % biom  R38
      '0.715 R5P + 3.624 AcCoA + 1.191 G6P + 0.501 E4P + 1.205 G3P + 1.002 PEP + 1.197 Pyr + 2.039 OAA+  1.233 AKG + 0.133 GAP -> Biomass + 0.683 Dummy_Fum.ex + 1.017 Dummy_CO2.ex';... % biom  R38
...
... % dummy reactions
... % ---------------------------------------------------------------------
      'OAA (abcd) + Dummy_Fum.ex -> Fum (abcd)';...      % R39
      'CO20.ex (a) -> CO2_bm (a)';...                    % R40
      'CO21.ex (a) -> CO2_bm (a)';...                    % R41
      'CO2_bm (a) + Dummy_CO2.ex -> CO2 (a)';...         % R40
});

% tracers
t = tracer({...
'13C_CO2: CO2.ex @ 1';...
'12C_CO2: GO2.ex @ ';...
'CO20: CO20.ex @ ';...
'CO21: CO21.ex @ 1';...
});
t.frac = [0.5,0.5,1,1];

% MS data: minimal
d_min = msdata({...
      'G3P3: G3P @ 1 2 3';...
});

% MS data: synthetic
d_sim = msdata({...
      'G3P2: G3P @ 2 3';...
      'G3P3: G3P @ 1 2 3';...
      'DHAP3: DHAP @ 1 2 3';...
      'PEP3: PEP @ 1 2 3';...
      'Fum4: Fum @ 1 2 3 4';...
      'R5P5: R5P @ 1 2 3 4 5';...
      'RuBP5: RuBP @ 1 2 3 4 5';...
      'Cit5: Cit @ 1 2 3 4 5';...
      'Cit6: Cit @ 1 2 3 4 5 6';...
      'S7P7: S7P @ 1 2 3 4 5 6 7';
});
d_real = msdata({...
      'G3P2: G3P @ 2 3';...
      'G3P3: G3P @ 1 2 3';...
      'DHAP3: DHAP @ 1 2 3';...
      'PEP3: PEP @ 1 2 3';...
      'Fum4: Fum @ 1 2 3 4';...
      'R5P5: R5P @ 1 2 3 4 5';...
      'RuBP5: RuBP @ 1 2 3 4 5';...
      'Cit5: Cit @ 1 2 3 4 5';...
      'Cit6: Cit @ 1 2 3 4 5 6';...
      'S7P7: S7P @ 1 2 3 4 5 6 7';...
      % additional measurements
      'F6P6: F6P @ 1 2 3 4 5 6';...
      'G6P6: G6P @ 1 2 3 4 5 6';...
      'GAP3: GAP @ 1 2 3';...
      'Ru5P5: Ru5P @ 1 2 3 4 5';...
      'Suc4: Suc @ 1 2 3 4';...
      'Mal3: Mal @ 2 3 4';...
      'Ga2: Ga @ 1 2';...
      'Ga3: Ga @ 1 2 3';
});

d.idvs = idv; % create devault IDVs, which will be replaced by simulated IDVs later.

% set up model
x = experiment(t);
x.data_ms = d_min;
%x.data_ms = d_sim;
%x.data_ms = d_real;
m = model(r);
m.expts = x;

% Take care of symmetrical metabolites
m.mets{'Suc'}.sym = list('rotate180',atommap('1:4 2:3 3:2 4:1'));
m.mets{'Fum'}.sym = list('rotate180',atommap('1:4 2:3 3:2 4:1'));

%% 

% simulate MS measurements
% INST
m.options.int_tspan = [20.0,40.0,60.0,90.0,130.0,250.0,480.0,610.0];   % time points in seconds from "real" data for comparison
m.options.sim_tunit = 'h';   % seconds are unit of time
m.options.fit_reinit = true;
m.options.sim_ss = false;
m.options.sim_sens = false;

m.options.int_reltol = 1e-9;

m.rates.flx.val = [...
      0.130835102594275,...     % pgi_f
      3.71668741215529,...      % pgi_b
      3.29363256816686,...      % g6pdh
      4.54398740647375,...      % pfk_f
      12.7881799451705,...      % pfk_b
      0.119341336137933,...     % fba_f
      8.36353387483463,...      % fba_b
      10.2934261293198,...      % tpi_f
      21.6700571677651,...      % tpi_b
      3.23124065939782,...      % gapdh_f
      29.9413539355077,...      % gapdh_b
      30.3945560336893,...      % gpm_f
      27.9245063654956,...      % gpm_b
      3.0886696473401,...       % eno_f
      0.618619979146446,...     % eno_b
      2.36091347157904,...      % pk_f
      1.18835443501743,...      % pk_b
      2.04616472470356,...      % rpe_f
      9.8369434535881,...       % rpe_b
      2.66961147708459,...      % rpi_f
      6.32810901926316,...      % rpi_b
      14.7429088392299,...      % prk
      14.7029088392299,...      % rbc1
      2.74832196577336,...      % tkt1_f
      10.5391006946579,...      % tkt1_b
      16.1553046187192,...      % tkt2_f
      12.1984534195,...         % tkt2_b
      4.90552152071357,...      % tkt3_f
      1.07159399104816,...      % tkt3_b
      214.281997662109,...      % tal1_f
      213.580508632192,...      % tal1_b
      30.019590124834,...       % tal2_f
      30.7210791547507,...      % tal2_b
      3.13243849974859,...      % sba
      3.13243849974859,...      % sbp
      1.37886715541236,...      % pdh
      0.406109731495618,...     % cs
      0.431992601739519,...     % can_f
      0.0258828702438167,...    % can_b
      0.678462140095176,...     % icd_f
      0.375937420415126,...     % icd_b
      2.85636472543462,...      % sdh_f
      2.752779713619,...        % sdh_b
      2.01148170584367,...      % fum_f
      1.74031811856861,...      % fum_b
      15.0590625784128,...      % mdh_f
      15.204313979322,...       % mdh_b
      0.103585011815645,...     % icl
      7.85043971907002,...      % ms_f
      7.76685470725435,...      % ms_b
      0.50000000000004,...      % me
      1.05164329252468,...      % ppc
      0.0399999999999112,...    % rbc2
      0.040000000000021,...     % pgp
      0.040000000000021,...     % gld
      0.0299999999999934,...    % gt
      4.30450476016747,...      % glyk_f
      4.27450476016738,...      % glyk_b
      10.0000000000001,...      % co2in
      0.245356625855699,...     % biom
      0.167578575459,...        % R39 (bm Fum)
      0.124763844,...           % R40 (bm CO20)
      0.124763844,...           % R41 (bm CO21)
      0.249527688,...           % R42 (bm CO2)
      
];
m.rates.flx.val = mod2stoich(m)'; % make sure the fluxes are feasible

m.states.lb = 1e-6;
m.states.ub = 1;
m.states.val = [...
      10,...       % AKG
      10,...       % AcCoA
      10,...       % Biomass
      10,...       % CO2
      10,...       % CO20
      10,...       % CO21
      10,...       % CO2_bm
      10,...       % CO2.ext ??? TODO
      10,...       % Cit
      10,...       % DHAP
      10,...       % Dummy_CO2
      10,...       % Dummy_Fum
      10,...       % E4P
      10,...       % F6P
      10,...       % FBP
      10,...       % Fum
      10,...       % G2P
      10,...       % G3P
      10,...       % G6P
      10,...       % GAP
      10,...       % Ga
      10,...       % Gc
      10,...       % Glx
      10,...       % ICit
      10,...       % Mal
      10,...       % OAA
      10,...       % PEP
      10,...       % PG
      10,...       % Pyr
      10,...       % R5P
      10,...       % Ru5P
      10,...       % RuBP
      10,...       % S7P
      10,...       % SBP
      10,...       % Suc
      10,...       % TA
      10,...       % TK
      10,...       % X5P
];

% simulate measurements
fprintf("Running pre-evaluation\n");
s = evalc('simulate(m)');
s = evalc('simulate(m)');

%maxNumCompThreads(1);
REPS=100;

fprintf("Doing %d evaluations\n", REPS);

time = zeros(REPS,1);
for i = 1:REPS
    [s, out] = evalc('simulate(m)');
    a = strsplit(s, "\n");
    b = strsplit(a{4},{':', ' '});
    time(i) = str2double(b{3});
end
disp("======================")
fprintf("runtime %f ms +- %f\n\n", 1000*mean(time), 1000*std(time));

fprintf("%f, ", 1000*time)
fprintf("\n")

