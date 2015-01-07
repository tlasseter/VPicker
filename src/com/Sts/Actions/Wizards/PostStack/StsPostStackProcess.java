package com.Sts.Actions.Wizards.PostStack;

import com.Sts.Actions.Wizards.PostStack2d.*;
import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.Vsp.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;

abstract public class StsPostStackProcess extends StsWizardStep
{
	public StsProgressPanel panel;
    protected StsHeaderPanel header;
	protected StsSeismicBoundingBox[] segyVolumes;
	protected StsSEGYFormatRec[] selectedAttributes;

	protected StsSeismicWizard wizard;
    protected String processString;
    protected int nSuccessfulVolumes = 0;

    abstract public void run();

    public int getNSuccessfulVolumes() { return nSuccessfulVolumes;}

    public StsPostStackProcess(StsSeismicWizard wizard)
	{
		super(wizard);
		this.wizard = wizard;
        processString = "Post-stack 3d";
	}

    public StsPostStackProcess(StsPostStack2dWizard wizard)
	{
		super(wizard);
		this.wizard = wizard;
        processString = "Post-stack 2d";
    }

    public StsPostStackProcess(StsSegyVspWizard wizard)
	{
		super(wizard);
		this.wizard = wizard;
        processString = "VSP volumes";
    }

    public boolean start()
	{
		run();
		return true;
	}

	public void constructPanel()
	{
		segyVolumes = wizard.getSegyVolumesToProcess();
		if (segyVolumes == null) return;

        panel = StsProgressPanel.constructorWithCancelButton(20, 50);
		header = new StsHeaderPanel();
		setPanels(panel, header);
        header.setTitle("Processing " + processString);
		header.setSubtitle("Process " + processString);
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#ProcessSeismic");
		header.setInfoText(wizardDialog, "(1) Once complete, press the Finish Button to dismiss the screen\n" +
								 "     ***** Proceed to the poststack seismic load workflow step to load and view the data *****\n" +
								 "     ***** Cancel will dismis the screen and stop the processing after the current line is complete *****\n");
        wizard.dialog.setTitle("Process " + processString);
	}

    protected void setCancelStatus(StsSeismicBoundingBox volume)
    {
        panel.progressBar.setDescriptionAndLevel("Cancelled", StsProgressBar.WARNING);
        if(volume == null)
            panel.appendLine("User has cancelled processing.");
        else
            panel.appendLine("User has cancelled processing for volume " + volume.getName());
    }

	public boolean end()
	{
		return true;
	}
}
