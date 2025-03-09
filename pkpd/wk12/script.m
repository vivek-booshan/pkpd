% less than 5% of food converted to disach by salivary amylase
% [St] / [S0] = 0.95
% t = 30 seconds of chewing
% -ln(0.95/30) = 0.00171 = kmouth
p.kmouth = 0.00171;
% 2 hours to get to the stomach and 50% dumped to intestine
% ln(0.5/120) = 0.0058 min^-1 --> 9.67e-5 s^-1
p.kFmouth_int = 9.67e-5;
p.kDmouth_int = 9.67e-5;
p.kFGint = 0.00064; % 90% conversion in 1 hr
p.kFMint = (4/6)*p.kFGint; % about 60% of food turn to gluc and rest is other
p.kpoop = 2e-7; % 12-24 hrs to poop and most sugars absorbed so basically zero
p.kDGint = p.kFGint; % same rate of conversion as food to glucose basically
p.kG_Gliver = 0.0033; % 5-15 minutes
p.kDMint = p.kFMint;
p.kGnG = 0.000142; % 40% of glucose converted in 1 hr, assuming no storage capacity
p.knGG = 1/p.kGnG; % doesn't really go back without glucagon signal
p.kG_Gblood = 0.0000115; % 100 g/day released
p.kGblood_Gtissue = 0.00001; 
p.kGtissue_Gblood = 1/0.00001;
p.kcl = 0.002; % renal threshold 180 mmol/dL
p.kcl2 = 0.0001; % muscle and fat and other stuff moving and consuming it

