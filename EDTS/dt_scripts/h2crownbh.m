function [crownbh]=h2crownbh(height,ipft)

b1Cl(1)     = 0.99;
b1Cl(2:4)   = 0.3106775;
b1Cl(5)     = 0.99;
b1Cl(6:11)  = 0.3106775;
b1Cl(12:16) = 0.99;
b1Cl(17)    = 0.3106775;

b2Cl(1)     = 1.0;
b2Cl(2:4)   = 1.098;
b2Cl(5)     = 1.0;
b2Cl(6:11)  = 1.098;
b2Cl(12:16) = 1.0;
b2Cl(17)    = 1.098;


crown_length = b1Cl(ipft) * height.^b2Cl(ipft);

crownbh    = max([0.00 height-crown_length]);


