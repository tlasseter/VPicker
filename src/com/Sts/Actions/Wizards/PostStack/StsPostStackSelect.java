package com.Sts.Actions.Wizards.PostStack;

import com.Sts.Actions.Wizards.PreStack3d.*;
import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.Vsp.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsPostStackSelect extends StsWizardStep
{
    public StsPostStackSelectPanel panel;
    public StsHeaderPanel header;
    boolean changeHeaders = false;

    public StsPostStackSelect(StsPostStackWizard wizard)
    {
        super(wizard);
        initializePanels();
        header.setTitle("Post-stack SegY Definition");
        header.setSubtitle("Select Post-Stack File(s)");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#ProcessSeismic");
        wizard.dialog.setTitle("Select Poststack File(s)");
        setInfoText();
    }

    public StsPostStackSelect(StsPreStackWizard wizard)
    {
        super(wizard);
        initializePanels();
        header.setTitle("Prestack SegY Definition");
        header.setSubtitle("Select Prestack File(s)");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#ProcessSeismic");
        wizard.dialog.setTitle("Select PresStack File(s)");
        setInfoText();
    }


    public StsPostStackSelect(StsSegyVspWizard wizard)
    {
        super(wizard);
        initializePanels();
        header.setTitle("VSP SegY Definition");
        header.setSubtitle("Select VSP Files(s)");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#ProcessSeismic");
        wizard.dialog.setTitle("Select VSP File(s)");
        setInfoText();
    }

    private void initializePanels()
    {
        panel = new StsPostStackSelectPanel(wizard, this);
        header = new StsHeaderPanel();
       setPanels(panel, header);
    }

    protected void setInfoText()
    {
        header.setInfoText(wizardDialog,"(1) Navigate to the directory containing the SegY Files using the Dir Button.\n" +
                                        "    All SegY Files in the selected directory will be placed in the left list.\n" +
                                        "      **** File names must have a .sgy, .SGY, .segy or .SEGY suffix *****\n" +
                                        "(2) Select files from the left list and place them in the right using the controls between the lists.\n" +
                                        "      **** All selected files will be scanned and placed in the table at bottom of screen with stats****\n" +
                                        "      **** from a random scan for review to ensure files are read correctly and disk space is adequate. ****\n" +
                                        "      **** If scan fails, adjust parameters on this and next screen before proceeding w/ processing ****\n" +
                                        "(3) Select a pre-defined Segy Format template if one exists for the selected files\n" +
                                        "      **** Format templates are saved from previous executions of this wizard.\n" +
                                        "(4) Set the scan percentage (# of samples of the SegY Files scanned) for verifying correct file reading.\n" +
                                        "      **** Random scanning of samples is used to determine if file is being read correctly.\n" +
                                        "(5) Once all selected file formats are correctly defined press the Next>> Button\n" +
                                        "      **** Dont worry if value in table are incorrect, next screen allows trace header mapping.****");
    }

    public StsPostStackSelectPanel getPanel()
    {
        return panel;
    }

    public boolean start()
    {
		panel.initialize();
        return true;
    }

    public boolean end()
    {
        try
        {
            if(panel.getVolumeName() == null)
            {
                new StsMessage(wizard.frame, StsMessage.ERROR, "No volumes selected or name specified: select or cancel.");
                return false;
            }
			StsSeismicWizard segyVolumeWizard = ((StsSeismicWizard) wizard);
            segyVolumeWizard.setVolumeName(panel.getVolumeName());
            segyVolumeWizard.setVolumesSegyFormat();
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsPostStackSelect.end() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    public void changeHeaders()
    {
        changeHeaders = true;
        wizard.next();
    }

    public boolean getToggleChangeHeaders()
    {
		changeHeaders = !changeHeaders; // toggle flag
		return !changeHeaders; // return original state
    }

	public void updatePanel()
	{
		panel.updatePanel();
	}

    public void clearFiles()
	{
		panel.clearFiles();
	}

	public void addAvailableVolume(StsSeismicBoundingBox volume)
	{
		panel.addAvailableVolume(volume);
	}

    public void removeVolume(StsSeismicBoundingBox volume)
    {
        panel.removeVolume(volume);
    }

    public void moveToAvailableList(StsSeismicBoundingBox volume)
    {
        panel.moveToAvailableList(volume);
    }
}
