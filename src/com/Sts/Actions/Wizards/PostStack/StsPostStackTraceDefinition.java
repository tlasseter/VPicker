package com.Sts.Actions.Wizards.PostStack;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.Interfaces.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2006</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */
abstract public class StsPostStackTraceDefinition extends StsWizardStep
{
	protected StsPostStackTraceDefinitionPanel panel;
	protected StsHeaderPanel header;
	protected StsSeismicWizard wizard = null;

    public StsPostStackTraceDefinition(StsWizard wizard)
	{
        super(wizard);
        this.wizard = (StsSeismicWizard)wizard;
        header = new StsHeaderPanel();
        header.setTitle(getTitle());
		header.setSubtitle(getSubtitle());
		//header.setLogo(getLogo());
		//header.setLink(getLink());
		header.setInfoText(wizardDialog, getInfoText().toString());
        panel = constructPanel();
        setPanels(panel, header);
    }

    abstract protected StsPostStackTraceDefinitionPanel constructPanel();

	public boolean start()
	{
		if(!panel.initialize()) return false;
//        panel.updatePanel();
        analyzeGrid();
        return true;
	}

	public void updatePanel()
	{
        StsToolkit.runLaterOnEventThread
        (
            new Runnable()
            {
                public void run()
                {
                    panel.updatePanel();
                }
            }
        );
    }
    /**
	 * getAllTraceRecords
	 *
	 * @param volume StsRotatedGridBoundingBox
	 * @return StsSEGYFormatRec[]
	 */
	public StsSEGYFormatRec[] getAllTraceRecords(StsSeismicBoundingBox volume)
	{
		return volume.getSegyFormat().getAllTraceRecords();
		//return wizard.getSegyFormat().getAllTraceRecords();
	}

	public StsSEGYFormatRec[] getRequiredTraceRecords(StsSeismicBoundingBox volume)
	{
            return volume.getSegyFormat().getRequiredTraceRecords();
	}

	abstract public StringBuffer getInfoText();

	public String getLink()
	{                      
		return "http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#ProcessSeismic";
	}

	public String getLogo()
	{
		return "S2SLogo.gif";
	}

	public int getNTraces(StsSeismicBoundingBox volume)
	{
		return volume.getTotalTraces();
	}

	public boolean getOverrideHeaderAllowed()
	{
		return false;
	}

	public String getSubtitle()
	{
		return "Trace Definition";
	}

	public String getTitle()
	{
		return "Post-stack Trace Header Details";
	}

	public byte[] getTraceHeaderBinary(StsSeismicBoundingBox volume, int nTrace)
	{
		return volume.getTraceHeaderBinary(nTrace);
	}

/**
	 * getTraceHeaderSize
	 *
	 * @param volume StsRotatedGridBoundingBox
	 * @return long
	 */
	public long getTraceHeaderSize(StsSeismicBoundingBox volume)
	{
		return volume.getSegyFormat().getTraceHeaderSize();
	}

	public boolean end()
	{
		return true;
	}

	public void initSetValues(StsSeismicBoundingBox volume)
	{
		volume.initializeSegyIO();
	}

	public int[] getTraceStartAndCount(StsSeismicBoundingBox volume, int traceNumber)
	{
		return new int[] {traceNumber, 24};
	}

	public boolean getOverrideStep()
	{
		return false;
	}

    /** This pattern allows the wizard to cancel the runnable; at the same time, the runnable has access
     *  to the wizard and can make make appropriate method calls.
     */
 /* Runnable must include the panelUpdate or it will be updated while runnable is still of on it's thread.  TJL 2/4/08
    public void analyzeGrid()
    {
        progressRunnable = StsPostStack3dAnalyzer.geometryAnalyzer(wizard, panel.progressPanel);
        StsToolkit.runRunnable(progressRunnable);
    }
 */
    public void analyzeGrid()
    {
         final StsPostStackAnalyzer analyzer = wizard.getAnalyzer(panel.progressPanel, panel.volumeStatusTablePanel);
         StsProgressRunnable progressRunnable = new StsProgressRunnable()
         {
              public void cancel()
             {
                 analyzer.cancelProcess();
             }
             public void run()
             {
                 if(!analyzer.analyzeGrid()) return;
                 updatePanel();
             }
         };
         StsToolkit.runRunnable(progressRunnable);
    }  
}