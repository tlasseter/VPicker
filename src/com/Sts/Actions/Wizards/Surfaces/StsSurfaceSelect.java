
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Surfaces;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.IO.*;

public class StsSurfaceSelect extends StsWizardStep
{
    StsSurfaceSelectPanel panel;
    StsHeaderPanel header;

    public StsSurfaceSelect(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSurfaceSelectPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Surface Selection");
        header.setSubtitle("Selecting Available Surfaces");
        header.setInfoText(wizardDialog,"(1) Select the surfaces that will be loaded.\n" +
                           "     ***** Supported formats can be found in the Users Guide ****\n" +
                           "(2) If other surfaces are desired, navigate to the directory containing the surfaces using the Dir button." +
                           "     ***** All surfaces in the selected directory will be placed in the left list. *****\n" +
                           "     ***** The default location is the project directory but any directory can be selected ****\n" +
                           "(3) Select the desired surfaces from the left list and place them in the right list using the provided controls between the lists.\n" +
                           "(4) Once all surface selections are complete, press the Next>> Button.");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Surfaces");

    }

//    public StsFile[] getSelectedFiles() {  return panel.getSelectedFiles(); }
    public StsAbstractFile[] getSelectedFiles() {  return panel.getSelectedFiles(); }
    public StsAbstractFile[] getViewableFiles() {  return panel.getViewableFiles(); }
    public String getSelectedDirectory() {  return panel.getCurrentDirectory(); }
    public boolean getReloadAscii() { return panel.getReloadAscii();  }

    public boolean start()
    {
        StsImportSeismicSurfaces.initializeColOrder();
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        return true;
    }
}

