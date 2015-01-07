package com.Sts.Actions.Wizards.PostStack;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Mar 21, 2008
 * Time: 9:25:51 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class StsPostStackAnalyzer
{
    protected StsSeismicWizard wizard;
    protected StsProgressPanel progressPanel;
    private StsTablePanelNew statusTablePanel;
    protected double scanPercent;
    protected StsSeismicBoundingBox[] selectedVolumes = null;
    protected StsSeismicBoundingBox[] totalVolumes = null;
    int nSelectedVolumes;
    int nTotalVolumes;
    protected boolean success = false;

    abstract public boolean analyzeGrid(StsProgressPanel progressPanel, boolean message, StsSeismicBoundingBox volume);

    public void initialize(StsSeismicWizard wizard, StsProgressPanel progressPanel, StsTablePanelNew statusTablePanel)
    {
        this.wizard = wizard;
        this.progressPanel = progressPanel;
        progressPanel.progressBar.canceled = false;
        this.statusTablePanel = statusTablePanel;
        this.selectedVolumes = wizard.getSelectedVolumes();
        nSelectedVolumes = selectedVolumes.length;
        this.totalVolumes = wizard.getSegyVolumes();
        nTotalVolumes = totalVolumes.length;
        scanPercent = wizard.getScanPercent();
        progressPanel.setLevel( StsProgressBar.INFO);
    }

    public boolean analyzeTraces(StsSeismicBoundingBox volume)
    {
        if(volume.analyzeTraces(scanPercent, progressPanel) != StsAnalyzeTraces.ANALYSIS_OK)
        {
            setFileStatus( StsSeismicBoundingBox.STATUS_TRACES_BAD, volume);
            setProgressErrorStatus("Failed to analyze traces for " + volume.getName() + ". ");
            return false;
        }
        setFileStatus(StsSeismicBoundingBox.STATUS_TRACES_OK, volume);
        setProgressStatus("Trace analysis " + volume.getName() + " OK.", StsProgressBar.INFO);
        return true;
    }

    public void cancelProcess()
    {
        progressPanel.cancel();
    }

    public void clearProcess()
    {
        progressPanel.clearProcess();
    }


    public synchronized boolean analyzeHeaders()
	{
        int nGoodVolumes = 0;
        try
		{
            if (selectedVolumes == null || selectedVolumes.length == 0)
			{
                setProgressErrorStatus("No volumes selected.");
				return false;
			}
            nGoodVolumes = 0;
            for (int n = 0; n < nTotalVolumes; n++)
                if(headerAnalysisOK(totalVolumes[n])) nGoodVolumes++;

            wizard.checkBlockButtons();
            progressPanel.initializeIntervals(nSelectedVolumes);
            for (int n = 0; n < nSelectedVolumes; n++)
			{
                progressPanel.setInterval(n);
                StsSeismicBoundingBox selectedVolume = selectedVolumes[n];
                if(headerAnalysisOK(totalVolumes[n])) continue;
                if(analyzeHeaders(selectedVolume)) nGoodVolumes++;
            }
            if(nGoodVolumes == 0)
                progressPanel.setLevel( StsProgressBar.ERROR);
            else if(nGoodVolumes < nTotalVolumes)
                progressPanel.setLevel( StsProgressBar.WARNING);

            progressPanel.finished();
            success = nGoodVolumes > 0;
            return success;
        }
		catch (Exception ex)
		{
			success = nGoodVolumes > 0;
			StsException.outputWarningException(this, "analyzeTraces() failed.", ex);
            setProgressErrorStatus("Exception thrown", ex.getMessage());
            return success;
        }
    }

    private boolean headerAnalysisOK(StsSeismicBoundingBox volume)
    {
        int status = volume.status;
        return status >= StsSeismicBoundingBox.STATUS_FILE_OK || status < StsSeismicBoundingBox.STATUS_FILE_BAD;
    }

    protected boolean analyzeHeaders(StsSeismicBoundingBox currentVolume)
    {
        try
		{
            if (isCanceled())
            {
                setProgressCancelStatus("User canceled", "User canceled.");
                return false;
            }

            currentVolume.readFileHeader();
            if (progressPanel.isCanceled())
            {
                setProgressCancelStatus("User canceled", "User canceled.");
                return false;
            }
            if (!currentVolume.analyzeBinaryHdr(progressPanel))
            {
                setFileStatus( StsSeismicBoundingBox.STATUS_HEADER_BAD, currentVolume);
                setProgressErrorStatus("Failed to analyze binary header for " + currentVolume.getName() + ". ");
                return false;
            }
            setFileStatus( StsSeismicBoundingBox.STATUS_HEADER_OK, currentVolume);
            if (!currentVolume.isAValidFileSize())
            {
                setFileStatus( StsSeismicBoundingBox.STATUS_FILE_BAD, currentVolume);
                setProgressErrorStatus("Invalid file size for " + currentVolume.getName() + ". ");
                return false;
            }
            setFileStatus( StsSeismicBoundingBox.STATUS_FILE_OK, currentVolume);
            if (progressPanel.isCanceled())
            {
                setProgressErrorStatus("User canceled.");
                return false;
            }
            return true;
        }
        catch (Exception ex)
		{
            setProgressErrorStatus("Exception thrown", ex.getMessage());
			StsException.outputWarningException(this, "analyzeHeaders", ex);
            return false;
        }
    }

    public synchronized boolean analyzeTraces()
	{
        StsSeismicBoundingBox currentVolume;
        int nVolumes;
        int nGoodVolumes = 0;
        try
		{
			if (selectedVolumes == null || selectedVolumes.length == 0)
			{
                setProgressErrorStatus("No volumes selected.");
				return false;
			}
            nGoodVolumes = 0;
            for (int n = 0; n < nTotalVolumes; n++)
                if(traceAnalysisOK(totalVolumes[n])) nGoodVolumes++;

            wizard.checkBlockButtons();
            progressPanel.initializeIntervals(nSelectedVolumes);
            for (int j = 0; j < nSelectedVolumes; j++)
			{
                progressPanel.setInterval(j);
                progressPanel.setDescriptionAndLevel("Analyzing traces for volume " + j, StsProgressBar.INFO);
                currentVolume = selectedVolumes[j];
                if(traceAnalysisOK(currentVolume))
                    continue;
                if(analyzeTraces(currentVolume)) nGoodVolumes++;
            }
            if(nGoodVolumes == 0)
                progressPanel.setLevel( StsProgressBar.ERROR);
            else if(nGoodVolumes < nTotalVolumes)
                progressPanel.setLevel( StsProgressBar.WARNING);
            progressPanel.finished();
            success = nGoodVolumes > 0;
            return success;
        }
		catch (Exception ex)
		{
			success = nGoodVolumes > 0;
			StsException.outputWarningException(this, "analyzeTraces", ex);
            setProgressErrorStatus("Exception thrown", ex.getMessage());
            return success;
        }
    }

    private boolean traceAnalysisOK(StsSeismicBoundingBox volume)
    {
        int status = volume.status;
        return status >= StsSeismicBoundingBox.STATUS_TRACES_OK || status < StsSeismicBoundingBox.STATUS_TRACES_BAD;

    }

    public boolean setStatus()
    {
        if (success)
            appendLine("All volumes analyzed successfully. ");
        else
            appendLine("Some or all volumes failed analysis.");
        wizard.checkUnblockButtons(success);
//      this repaint shouldn't be needed, but can't get the statusTablePanel to repaint properly without it
        wizard.dialog.repaint();
        return success;
    }

    protected boolean setGeometryStatus()
    {
        if (success)
            appendLine("All volumes analyzed successfully. ");
        else
            appendLine("Some or all volumes failed analysis.");
        wizard.checkUnblockButtons(success || wizard.getOverrideGeometry());
//      this repaint shouldn't be needed, but can't get the statusTablePanel to repaint properly without it
        wizard.dialog.repaint();
        return success;
    }

    public boolean analyzeGrid()
	{
		StsSeismicBoundingBox currentVolume;
        int nGoodVolumes = 0;
        try
		{
            success = true;
			if (selectedVolumes == null || selectedVolumes.length == 0)
			{
                setProgressErrorStatus("No volumes selected.");
                success = false;
                return false;
			}
            // TODO: Seems to be a problem here. We need to run the postStackAnalyzer on all volumes for which gridAnalysisOK.
            // TODO: instead it was being run only on selected volumes which is just one.
            // TODO: deleted two lines below and modified two lines below that. TJL 7/9/09
            nGoodVolumes = 0;
            // for (int n = 0; n < nTotalVolumes; n++)
            //     if(gridAnalysisOK(totalVolumes[n])) nGoodVolumes++;

            wizard.checkBlockButtons();
            // progressPanel.initialize(nSelectedVolumes);
            // for (int j = 0; j < nSelectedVolumes; j++)
            progressPanel.initialize(nTotalVolumes);
            for (int j = 0; j < nTotalVolumes; j++)
			{
                currentVolume = totalVolumes[j];
                progressPanel.setValue(j+1);
                if(gridAnalysisOK(currentVolume))
                    continue;
				if (analyzeGrid(progressPanel, true, currentVolume))
                    nGoodVolumes++;
            }
            if(nGoodVolumes == 0)
                progressPanel.setLevel( StsProgressBar.ERROR);
            else if(nGoodVolumes < nTotalVolumes)
                progressPanel.setLevel( StsProgressBar.WARNING);
            progressPanel.finished();
            success = nGoodVolumes > 0;
            return success;
        }
		catch (Exception ex)
		{
            success = nGoodVolumes > 0;
            setProgressErrorStatus("Exception thrown", ex.getMessage());
			StsException.outputWarningException(this, "analyzeTraces", ex);
            return success;
        }
        finally
        {
            setStatus();
        }
    }

    private boolean gridAnalysisOK(StsSeismicBoundingBox volume)
    {
        int status = volume.status;
        return status >= StsSeismicBoundingBox.STATUS_GEOMETRY_OK;
    }
    protected boolean isCanceled()
    {
        return progressPanel.isCanceled();
    }

    private void appendLine(String line)
    {
        progressPanel.appendLine(line);
    }

    protected void setProgressCancelStatus(String barMessage, String textMessage)
    {
        success = false;
        setProgressStatus( barMessage, textMessage, StsProgressBar.WARNING);
    }

    protected void setProgressErrorStatus(String barMessage, String textMessage)
    {
        success = false;
        progressPanel.finished();
        progressPanel.setDescriptionAndLevel(barMessage, StsProgressBar.ERROR);
        progressPanel.appendErrorLine(textMessage);
    }

    protected void setProgressErrorStatus(String textMessage)
    {
        success = false;
        progressPanel.finished();
        progressPanel.setDescriptionAndLevel(textMessage, StsProgressBar.ERROR);
        progressPanel.appendErrorLine(textMessage);
    }

    private void setProgressStatus(String barMessage, String textMessage, int level)
    {
        progressPanel.setDescriptionAndLevel(barMessage, level);
        progressPanel.appendLine(textMessage);
    }

    protected void setProgressStatus(String textMessage, int level)
    {
        progressPanel.setDescriptionAndLevel(textMessage, level);
        progressPanel.appendLine(textMessage);
    }

    protected void setFileStatus( int status, StsSeismicBoundingBox volume)
    {
        if( volume != null)
        {
            volume.status = status;
            statusTablePanel.setValueAt( StsSeismicBoundingBox.statusText[volume.status], volume, "statusString");
        }
    }
}
