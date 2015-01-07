
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Well;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

public class StsWellSelect extends StsWizardStep
{
    StsWellSelectPanel panel;
    StsHeaderPanel header, info;

    public StsWellSelect(StsWizard wizard, boolean reloadAscii)
    {
        super(wizard);
        panel = new StsWellSelectPanel(wizard, this, reloadAscii);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Well Selection");
        header.setSubtitle("Selecting Available Wells");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Well");
        header.setInfoText(wizardDialog,"(1) Navigate to the directory containing the wells using the Dir button. \n" +
                           " **** All wells in supported formats in the selected directory will be placed in the left list.\n" +
                           " ****     S2S Format (md, deltaX, deltaY) + logs + markers + time-depth\n" +
                           " ****     LAS Format (md, deltaX, deltaY) + logs\n" +
                           " ****     UT Format (md, drift, azimuth)\n" +
                           " ****     WLS Geographix formatted files\n" +
                           " **** If not geographix, any combination of formats can be used as long as they adhere to the S2S naming convention\n" +
                           " ****     <type>.<format>.<name>.<version>\n" +
                           " ****     <type> = well-dev, well-logs, well-ref, and well-td\n" +
                           " ****     <format> = txt, las, ut\n" +
                           " ****     <name> = User defined and must be same for all related files.\n" +
                           " ****     <version> = Optional integer to allow version control\n" +
                           "(3) Select the desired wells and place them in the right list using the provided controls\n " +
                           " **** Only one wls file can be selected and loaded at a time, any number of wells can be contained in the file.\n" +
                           "(4) Well file text can be viewed prior to loading by pressing the View File Button\n" +
                           "(5) Specify datum shift to apply to all read MD and Depth values. Must reload from ASCII to apply.\n" +
                           "(6) Once well selections are complete, press the Next>> Button.");
    }

//    public StsFile[] getSelectedFiles() {  return panel.getSelectedFiles(); }
    public String[] getSelectedWells() {  return panel.getSelectedWells(); }
    public String getSelectedDirectory() {  return panel.getSelectedDirectory(); }
//    public Object[] getSelectedWellObjects() {  return panel.getSelectedWellObjects(); }
//    public boolean getReloadASCII() { return panel.getReloadASCII(); }
//    public void setReloadASCII(boolean value) { panel.setReloadASCII(value); }

    public boolean start()
    {
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        return true;
    }
}

