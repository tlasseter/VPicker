package com.Sts.Actions.Wizards.Vsp;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;

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

public class StsPostStackVspAnalyzer extends StsPostStackAnalyzer
{
    static StsPostStackVspAnalyzer analyzer = null;

    static public StsPostStackVspAnalyzer getAnalyzer(StsSegyVspWizard wizard, StsProgressPanel progressPanel, StsTablePanelNew volumeStatusPanel)
    {
        if(analyzer != null)
        {
            // cancel any running analyzer; analyze routines are synchronized,
            // so new analyzer will wait until any old one is finished
            analyzer.cancelProcess();
        }
        else
            analyzer = new StsPostStackVspAnalyzer();
        analyzer.initialize(wizard, progressPanel, volumeStatusPanel);
        return analyzer;
    }

    public boolean analyzeGrid(StsProgressPanel progressPanel, boolean messages, StsSeismicBoundingBox volume)
    {
        StsSegyVsp currentVolume = (StsSegyVsp)volume;
        if(!currentVolume.analyzeGrid(progressPanel, messages))
        {
            setFileStatus( StsSeismicBoundingBox.STATUS_GRID_BAD, currentVolume);
            setProgressErrorStatus("Failed to analyze vsp geometry for " + currentVolume.getName() + ". ");
            return false;
        }
        else
        {
            setFileStatus( StsSeismicBoundingBox.STATUS_GEOMETRY_OK, currentVolume);
            setProgressStatus("PostStack2d " + currentVolume.getName() + " OK", StsProgressBar.INFO);
            return true;
        }
    }
}
