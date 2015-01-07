package com.Sts.Actions.Wizards.PostStack;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.Utilities.*;

abstract public class StsPostStackFileFormat extends StsWizardStep
{
	protected StsPostStackFileFormatPanel panel;
	protected StsHeaderPanel header;
    protected StsSeismicWizard wizard;

    public StsPostStackFileFormat(StsSeismicWizard wizard)
	{
		super(wizard);
        this.wizard = wizard;
        panel = new StsPostStackFileFormatPanel(wizard, this);
		header = new StsHeaderPanel();
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#ProcessSeismic");
		setPanels(panel, header);
        setTitlesAndInfo();
	}

    abstract protected void setTitlesAndInfo();

    public StsPostStackFileFormatPanel getPanel()
	{
		return panel;
	}

	public boolean start()
	{
        wizard.dialog.setTitle("Check and Set File Formats");
        panel.initialize();
//        reanalyzeHeaders();
        return true;
	}

	public boolean end()
	{
        return true;
	}

    public void reanalyzeHeaders()
    {

        Runnable runnable = new Runnable()
        {
             public void run()
             {
                StsPostStackAnalyzer analyzer = wizard.getAnalyzer(panel.progressPanel, panel.volumeStatusTablePanel);
                updateFileFormatPanel();
                analyzer.analyzeHeaders();
                analyzer.analyzeTraces();
                analyzer.setStatus();
                updatePanel();
             }
        };
        StsToolkit.runRunnable(runnable);
    }

    public void reanalyzeTraces()
    {

        Runnable runnable = new Runnable()
        {
             public void run()
             {
                StsPostStackAnalyzer analyzer = wizard.getAnalyzer(panel.progressPanel, panel.volumeStatusTablePanel);
                analyzer.analyzeTraces();
                analyzer.setStatus();
                updatePanel();
             }
        };
        StsToolkit.runRunnable(runnable);
    }

    public void updateFileFormatPanel()
    {
        panel.updateFileFormatPanel();
    }

	public void updatePanel()
	{
		panel.updatePanel();
	}
}