function dydt = bloodSugarODE(t, y, p)
    
    % F = food
    % D = disachharide
    % G = glucose
    % M = monosachharide (not glucose)
    % nG = glycogen

    dydt = zeros(11, 1);
    Fmouth = y(1);
    Dmouth = y(2);
    Fintestine = y(3);
    Dintestine = y(4);
    Gintestine = y(5);
    Mintestine = y(6);
    Gliver = y(7);
    nGliver = y(8);
    Gblood = y(9);
    Gtissue = y(10);
    
    % kmouth : food to disaccharide
    % kFmouth_int : food from mouth/stomach to intestine
    % kDmouth_int : disaccharide from mouth/stomach to intestine
    % kintestineFD : food to disaccharide
    % kpoop : poop
    % kDint_nGliver : disaccharide in intestine to (not Glucose) liver
    % kDGint : disaccharide to glucose in the intestine
    % kG_Gliver : glucose from intestine to liver
    % kDMint : disaccharide to monosaccharide within intestine
    % kGnG : glucose to not glucose
    % knGG : not glucose to glucose
    % kG_Gblood : liver glucose to blood glucose
    % kGblood_Gtissue : blood glucose to tissue glucose
    % kGtissue_Gblood : tissue glucose to blood glucose
    % kcl : kidney clearance
    % kcl2 : kidney clearance

    %%%%%%% MOUTH %%%%%%%%%
    
    dFmouthdt = -p.kmouth * 0.85 * Fmouth - p.kFmouth_int * 0.90 * Fmouth;
    dDmouthdt = p.kmouth * 0.85 * Fmouth - p.kDmouth_int * Dmouth;
    
    %%%%%%% INTESTINE %%%%%
    
    dFintestinedt = p.kFmouth_int * Fmouth -p.kFMint * Fintestine - p.kFGint * Fintestine;
    dDintestinedt = p.kDmouth_int * Dmouth - p.kDGint * Dintestine - p.kDMint * Dintestine; % nearly all becomes monosach
    dGintestinedt = p.kDGint * Dintestine - p.kG_Gliver * Gintestine; 
    dMintestinedt = p.kDMint * Dintestine;

    %%%%%%% LIVER %%%%%%%%%

    dGliverdt = p.kG_Gliver * Gintestine - p.kGnG * Gliver + p.nGG * nGliver - p.kG_Gblood * Gliver;
    dnGliverdt = p.kGnG * Gliver - p.nGG * nGliver;

    %%%%%%% BLOOD %%%%%%%%%

    dGblooddt = p.kG_Gblood * Gliver - p.kGblood_Gtissue * Gblood + p.kGtissue_Gblood * Gtissue - p.kcl * Gblood;

    %%%%%%% TISSUE %%%%%%%%

    dGtissuedt = p.kGblood_Gtissue * Gblood - p.kGtissue_Gblood * Gtissue - p.kcl2 * Gtissue;
    dydt = [dFmouthdt, dDmouthdt, dFintestinedt, dDintestinedt, dGintestinedt, dMintestinedt, dGliverdt, dnGliverdt, dGblooddt, dGtissuedt];

end