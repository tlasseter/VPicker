package com.Sts.Actions.Wizards.Well;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.UI.*;

public class StsOtherFiles extends StsWizardStep
{
    StsOtherFilesPanel panel;
    StsHeaderPanel header;

    public StsOtherFiles(StsWizard wizard)
    {
        super(wizard);
        panel = new StsOtherFilesPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
//        panel.setPreferredSize(new Dimension(500, 500));
        header.setTitle("Well Load");
        header.setSubtitle("Other Well File Definition");
        header.setInfoText(wizardDialog,"(1) \n" +
                           "(2) Press the Next>> Button.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Well");
    }

    public boolean start()
    {
        return true;
    }

    public boolean end()
    {
        if(panel.getSurveyFile() == null)
        {
            new StsMessage(wizard.frame, StsMessage.WARNING, "Survey file is required to position the well in 3D.");
            return false;
        }

        // Make sure the wells in the header file exist in the selected files.
        String[] wellNames = StsMultiWellImport.getWellNamesFromHeaderFile(((StsWellWizard)wizard).getFilename(StsWellWizard.HEADER),
                ((StsWellWizard)wizard).getNumberRowsToSkip(StsWellWizard.HEADER), ", ");

        // Examine Survey
        boolean inFile = true;
        String fileName = panel.getSurveyFile();
        for(int i=0; i<wellNames.length; i++)
        {
            inFile = StsMultiWellImport.isWellInFile(fileName,wellNames[i]);
            if(!inFile)
            {
                if(!StsYesNoDialog.questionValue(wizard.frame,
                    "Well " + wellNames[i] + " is not in selected survey file.\nContinue (Yes) or Cancel (No)?"))
                    return false;
            }
        }
        // Examine Time-Depth
        inFile = true;
        fileName = panel.getTdFile();
        if(fileName != null)
        {
            for(int i=0; i<wellNames.length; i++)
            {
                inFile = StsMultiWellImport.isWellInFile(fileName,wellNames[i]);
                if(!inFile)
                {
                    if(!StsYesNoDialog.questionValue(wizard.frame,
                        "Well " + wellNames[i] + " is not in selected time-depth file.\nContinue (Yes) or Cancel (No)?"))
                        return false;
                }
            }
        }
        // Examine Tops
        inFile = true;
        fileName = panel.getTopsFile();
        if(fileName != null)
        {
            for(int i=0; i<wellNames.length; i++)
            {
                inFile = StsMultiWellImport.isWellInFile(fileName,wellNames[i]);
                if(!inFile)
                {
                    if(!StsYesNoDialog.questionValue(wizard.frame,
                        "Well " + wellNames[i] + " is not in selected tops file.\nContinue (Yes) or Cancel (No)?"))
                        return false;
                }
            }
        }
        // Examine Log Files
        inFile = false;
        String[] logFilenames = panel.getLogFiles();
        if(logFilenames == null)
            return true;

        for(int i=0; i<wellNames.length; i++)
        {
            for(int j=0; j<logFilenames.length; j++)
            {
                inFile = StsMultiWellImport.isWellInFile(logFilenames[j],wellNames[i]);
                if(inFile)
                    break;
            }
            if(!inFile)
            {
                 if(!StsYesNoDialog.questionValue(wizard.frame,
                        "Well " + wellNames[i] + " is not in any of the selected log files.\nContinue (Yes) or Cancel (No)?"))
                 return false;
            }
        }
        return true;
    }

}
