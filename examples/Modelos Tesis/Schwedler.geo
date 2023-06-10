// Gmsh project created on Tue Mar 21 20:12:25 2023
//+
SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 22.9, 0, 2*Pi};
//+
Point(2) = {0, 0, 0, 1.0};
//+
Point(3) = {0, 22.9, 0, 1.0};
//+
Line(2) = {2, 3};
//+//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{2}; Point{3}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{3}; Point{3}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{4}; Point{7}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{5}; Point{9}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{6}; Point{11}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{7}; Point{13}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{8}; Point{15}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{9}; Point{17}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{10}; Point{19}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{11}; Point{21}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{12}; Point{23}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Point{25}; Curve{13}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{14}; Point{27}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Point{29}; Curve{15}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Point{31}; Curve{16}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Point{33}; Curve{17}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Point{35}; Curve{18}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Point{37}; Curve{19}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Point{39}; Curve{20}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Point{1}; Curve{21}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Point{43}; Curve{22}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Point{45}; Curve{23}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Point{47}; Curve{24}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Point{49}; Curve{25}; }
}
//+
Circle(27) = {0, 0, 1.79, 21.37, 0, 2*Pi};
//+
Circle(28) = {0, 0, 3.26, 16.41, 0, 2*Pi};
//+
Circle(29) = {0, 0, 4.27, 8.78, 0, 2*Pi};
//+
Point(55) = {0, 0, 4.58, 1.0};
//+
Line(30) = {1, 52};
//+
Line(31) = {52, 53};
//+
Line(32) = {53, 54};
//+
Line(33) = {54, 55};
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{30}; Point{52}; Curve{31}; Point{53}; Curve{32}; Point{54}; Curve{33}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{34}; Point{62}; Curve{35}; Point{63}; Curve{36}; Point{64}; Curve{37}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{38}; Point{72}; Curve{39}; Point{73}; Curve{40}; Point{74}; Curve{41}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{42}; Point{82}; Curve{43}; Point{83}; Curve{44}; Point{84}; Curve{45}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{46}; Point{92}; Curve{47}; Point{93}; Curve{48}; Point{94}; Curve{49}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{50}; Point{102}; Curve{51}; Point{103}; Curve{52}; Point{104}; Curve{53}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{54}; Point{112}; Curve{55}; Point{113}; Curve{56}; Point{114}; Curve{57}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{58}; Point{122}; Curve{59}; Point{123}; Curve{60}; Point{124}; Curve{61}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{62}; Point{132}; Curve{63}; Point{133}; Curve{64}; Point{134}; Curve{65}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{66}; Point{142}; Curve{67}; Point{143}; Curve{68}; Curve{22}; Point{144}; Curve{69}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{8}; Curve{70}; Point{154}; Curve{71}; Point{155}; Curve{72}; Point{156}; Curve{74}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{76}; Point{170}; Curve{77}; Point{171}; Curve{78}; Point{172}; Curve{79}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{80}; Point{180}; Curve{81}; Point{181}; Curve{82}; Point{182}; Curve{83}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{11}; Curve{84}; Point{190}; Curve{85}; Point{191}; Curve{86}; Point{192}; Curve{87}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{89}; Point{204}; Curve{90}; Point{205}; Curve{91}; Point{206}; Curve{92}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{93}; Point{214}; Curve{94}; Point{215}; Curve{95}; Point{216}; Curve{96}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{97}; Point{224}; Curve{98}; Point{225}; Curve{99}; Point{226}; Curve{100}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{101}; Point{234}; Curve{102}; Point{235}; Curve{103}; Point{236}; Curve{104}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{105}; Point{244}; Curve{106}; Point{245}; Curve{107}; Point{246}; Curve{108}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{109}; Point{254}; Curve{110}; Point{255}; Curve{111}; Point{256}; Curve{112}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{113}; Point{264}; Curve{114}; Point{265}; Curve{115}; Point{266}; Curve{116}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{117}; Point{274}; Curve{118}; Point{275}; Curve{119}; Point{276}; Curve{120}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 1}, 2*Pi/24} {
  Duplicata { Curve{121}; Point{284}; Curve{122}; Point{285}; Curve{123}; Point{286}; Curve{124}; }
}
//+
Circle(129) = {0, 0, 0, 22.9, 0, 2*Pi};
//+
Recursive Delete {
  Point{2}; 
}
//+
Recursive Delete {
  Point{2}; 
}
//+
Recursive Delete {
  Point{2}; 
}
//+
Recursive Delete {
  Curve{18}; Curve{17}; Curve{16}; Curve{15}; Curve{14}; Curve{13}; Curve{12}; Curve{11}; Curve{10}; Curve{9}; Curve{8}; Curve{7}; Curve{6}; Curve{5}; Curve{4}; Curve{3}; Curve{2}; Curve{25}; Curve{24}; Curve{23}; Curve{22}; Curve{21}; Curve{20}; Curve{19}; 
}
//+
Recursive Delete {
  Curve{73}; Curve{26}; 
}
//+
Recursive Delete {
  Curve{88}; Curve{75}; 
}
//+
Line(130) = {101, 92};
//+
Line(131) = {102, 93};
//+
Line(132) = {103, 94};
//+
Line(133) = {91, 82};
//+
Line(134) = {92, 83};
//+
Line(135) = {93, 84};
//+
Line(136) = {81, 72};
//+
Line(137) = {82, 73};
//+
Line(138) = {83, 74};
//+
Line(139) = {71, 62};
//+
Line(140) = {72, 63};
//+
Line(141) = {73, 64};
//+
Line(142) = {61, 52};
//+
Line(143) = {62, 53};
//+
Line(144) = {63, 54};
//+
Line(145) = {1, 294};
//+
Line(146) = {52, 295};
//+
Line(147) = {53, 296};
//+
Line(148) = {293, 284};
//+
Line(149) = {294, 285};
//+
Line(150) = {295, 286};
//+
Line(151) = {283, 274};
//+
Line(152) = {284, 275};
//+
Line(153) = {285, 276};
//+
Line(154) = {273, 264};
//+
Line(155) = {274, 265};
//+
Line(156) = {275, 266};
//+
Line(157) = {263, 254};
//+
Line(158) = {264, 255};
//+
Line(159) = {265, 256};
//+
Line(160) = {253, 244};
//+
Line(161) = {254, 245};
//+
Line(162) = {255, 246};
//+
Line(163) = {243, 234};
//+
Line(164) = {244, 235};
//+
Line(165) = {245, 236};
//+
Line(166) = {233, 224};
//+
Line(167) = {234, 225};
//+
Line(168) = {235, 226};
//+
Line(169) = {225, 216};
//+
Line(170) = {224, 215};
//+
Line(171) = {223, 214};
//+
Line(172) = {213, 204};
//+
Line(173) = {214, 205};
//+
Line(174) = {215, 206};
//+
Line(175) = {203, 190};
//+
Line(176) = {204, 191};
//+
Line(177) = {205, 192};
//+
Line(178) = {189, 180};
//+
Line(179) = {190, 181};
//+
Line(180) = {191, 182};
//+
Line(181) = {179, 170};
//+
Line(182) = {180, 171};
//+
Line(183) = {181, 172};
//+
Line(184) = {171, 156};
//+
Line(185) = {170, 155};
//+
Line(186) = {169, 154};
//+
Line(187) = {153, 142};
//+
Line(188) = {154, 143};
//+
Line(189) = {155, 144};
//+
Line(190) = {143, 134};
//+
Line(191) = {142, 133};
//+
Line(192) = {141, 132};
//+
Line(193) = {131, 122};
//+
Line(194) = {132, 123};
//+
Line(195) = {133, 124};
//+
Line(196) = {121, 112};
//+
Line(197) = {122, 113};
//+
Line(198) = {123, 114};
//+
Line(199) = {111, 102};
//+
Line(200) = {112, 103};
//+
Line(201) = {113, 104};
//+//+
Recursive Delete {
  Curve{27}; 
}
//+
Recursive Delete {
  Curve{28}; 
}
//+
Recursive Delete {
  Curve{29}; 
}
//+
Recursive Delete {
  Curve{1}; 
}
//+
Recursive Delete {
  Curve{129}; 
}
//+
Line(202) = {203, 213};
//+
Line(203) = {213, 223};
//+
Line(204) = {223, 233};
//+
Line(205) = {233, 243};
//+
Line(206) = {243, 253};
//+
Line(207) = {253, 263};
//+
Line(208) = {263, 273};
//+
Line(209) = {273, 283};
//+
Line(210) = {283, 293};
//+
Line(211) = {293, 1};
//+
Line(212) = {1, 61};
//+
Line(213) = {61, 71};
//+
Line(214) = {71, 81};
//+
Line(215) = {81, 91};
//+
Line(216) = {91, 101};
//+
Line(217) = {101, 111};
//+
Line(218) = {111, 121};
//+
Line(219) = {121, 131};
//+
Line(220) = {131, 141};
//+
Line(221) = {141, 153};
//+
Line(222) = {153, 169};
//+
Line(223) = {169, 179};
//+
Line(224) = {179, 189};
//+
Line(225) = {189, 203};
//+
Line(226) = {122, 132};
//+
Line(227) = {132, 142};
//+
Line(228) = {142, 154};
//+
Line(229) = {154, 170};
//+
Line(230) = {170, 180};
//+
Line(231) = {180, 190};
//+
Line(232) = {190, 204};
//+
Line(233) = {204, 214};
//+
Line(234) = {214, 224};
//+
Line(235) = {224, 234};
//+
Line(236) = {234, 244};
//+
Line(237) = {244, 254};
//+
Line(238) = {254, 264};
//+
Line(239) = {264, 274};
//+
Line(240) = {274, 284};
//+
Line(241) = {284, 294};
//+
Line(242) = {294, 52};
//+
Line(243) = {52, 62};
//+
Line(244) = {62, 72};
//+
Line(245) = {72, 82};
//+
Line(246) = {82, 92};
//+
Line(247) = {92, 102};
//+
Line(248) = {102, 112};
//+
Line(249) = {112, 122};
//+
Line(250) = {63, 73};
//+
Line(251) = {73, 83};
//+
Line(252) = {83, 93};
//+
Line(253) = {93, 103};
//+
Line(254) = {103, 113};
//+
Line(255) = {113, 123};
//+
Line(256) = {123, 133};
//+
Line(257) = {133, 143};
//+
Line(258) = {143, 155};
//+
Line(259) = {155, 171};
//+
Line(260) = {171, 181};
//+
Line(261) = {181, 191};
//+
Line(262) = {191, 205};
//+
Line(263) = {205, 215};
//+
Line(264) = {215, 225};
//+
Line(265) = {225, 235};
//+
Line(266) = {235, 245};
//+
Line(267) = {245, 255};
//+
Line(268) = {255, 265};
//+
Line(269) = {265, 275};
//+
Line(270) = {275, 285};
//+
Line(271) = {285, 295};
//+
Line(272) = {295, 53};
//+
Line(273) = {53, 63};
//+
Line(274) = {266, 276};
//+
Line(275) = {276, 286};
//+
Line(276) = {286, 296};
//+
Line(277) = {296, 54};
//+
Line(278) = {54, 64};
//+
Line(279) = {64, 74};
//+
Line(280) = {74, 84};
//+
Line(281) = {84, 94};
//+
Line(282) = {94, 104};
//+
Line(283) = {104, 114};
//+
Line(284) = {114, 124};
//+
Line(285) = {124, 134};
//+
Line(286) = {134, 144};
//+
Line(287) = {144, 156};
//+
Line(288) = {156, 172};
//+
Line(289) = {172, 182};
//+
Line(290) = {182, 192};
//+
Line(291) = {192, 206};
//+
Line(292) = {206, 216};
//+
Line(293) = {216, 226};
//+
Line(294) = {226, 236};
//+
Line(295) = {236, 246};
//+
Line(296) = {246, 256};
//+
Line(297) = {256, 266};
//+
Recursive Delete {
  Curve{87}; 
}
//+
Recursive Delete {
  Curve{92}; 
}
//+
Recursive Delete {
  Curve{96}; 
}
//+
Recursive Delete {
  Curve{100}; 
}
//+
Recursive Delete {
  Curve{104}; 
}
//+
Recursive Delete {
  Curve{108}; 
}
//+
Recursive Delete {
  Curve{112}; 
}
//+
Recursive Delete {
  Curve{116}; 
}
//+
Recursive Delete {
  Curve{120}; 
}
//+
Recursive Delete {
  Curve{124}; 
}
//+
Recursive Delete {
  Curve{128}; 
}
//+
Recursive Delete {
  Curve{33}; 
}
//+
Recursive Delete {
  Curve{37}; 
}
//+
Recursive Delete {
  Curve{41}; 
}
//+
Recursive Delete {
  Curve{45}; 
}
//+
Recursive Delete {
  Curve{49}; 
}
//+
Recursive Delete {
  Curve{53}; 
}
//+
Recursive Delete {
  Curve{57}; 
}
//+
Recursive Delete {
  Curve{61}; 
}
//+
Recursive Delete {
  Curve{65}; 
}
//+
Recursive Delete {
  Curve{69}; 
}
//+
Recursive Delete {
  Curve{74}; 
}
//+
Recursive Delete {
  Curve{79}; 
}
//+
Point(297) = {0, 0, 4.58, 1.0};
//+
Line(298) = {182, 297};
//+
Line(299) = {172, 297};
//+
Line(300) = {156, 297};
//+
Line(301) = {144, 297};
//+
Line(302) = {134, 297};
//+
Line(303) = {124, 297};
//+
Line(304) = {114, 297};
//+
Line(305) = {104, 297};
//+
Line(306) = {94, 297};
//+
Line(307) = {84, 297};
//+
Line(308) = {74, 297};
//+
Line(309) = {64, 297};
//+
Line(310) = {54, 297};
//+
Line(311) = {296, 297};
//+
Line(312) = {286, 297};
//+
Line(313) = {276, 297};
//+
Line(314) = {266, 297};
//+
Line(315) = {256, 297};
//+
Line(316) = {246, 297};
//+
Line(317) = {236, 297};
//+
Line(318) = {226, 297};
//+
Line(319) = {216, 183};
//+
Line(320) = {206, 183};
//+
Recursive Delete {
  Curve{320}; 
}
//+
Recursive Delete {
  Curve{319}; 
}
//+
Recursive Delete {
  Curve{318}; 
}
//+
Recursive Delete {
  Curve{83}; 
}
//+
Recursive Delete {
  Curve{298}; 
}
//+
Line(318) = {226, 297};
//+
Line(319) = {216, 297};
//+
Line(320) = {206, 297};
//+
Line(321) = {192, 297};
//+
Line(322) = {182, 297};

Physical Point("00_01_01_00") = {1,61,71,81,91,101,111,121,131,141,153,169,179,189,203,213,223,233,243,253,263,273,283,293};

Physical Point ("00_01_02_00") = {297};

Physical Line ("01_02_00_00") = {30:322};
