package com.Sts.Actions.Wizards.PostStack3d;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.Actions.Wizards.Seismic.*;
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

public class StsPostStack3dAnalyzer extends StsPostStackAnalyzer
{
    private static StsPostStack3dAnalyzer analyzer;

    static public StsPostStack3dAnalyzer getAnalyzer(StsSeismicWizard wizard, StsProgressPanel progressPanel, StsTablePanelNew volumeStatusPanel)
    {
        if(analyzer != null)
        {
            // cancel any running analyzer; analyze routines are synchronized,
            // so new analyzer will wait until any old one is finished
            analyzer.clearProcess();
        }
        else
            analyzer = new StsPostStack3dAnalyzer();
        analyzer.initialize(wizard, progressPanel, volumeStatusPanel);
        return analyzer;
    }

    public boolean analyzeGrid(StsProgressPanel progressPanel, boolean messages, StsSeismicBoundingBox volume)
    {
        StsSegyVolume currentVolume = (StsSegyVolume)volume;
        if(!currentVolume.analyzeGrid(progressPanel, messages))
        {
            setFileStatus( StsSeismicBoundingBox.STATUS_GRID_BAD, currentVolume);
            setProgressErrorStatus("Failed to analyze iline/xline geometry for " + currentVolume.getName() + ". ");
            return false;
        }
        setFileStatus( StsSeismicBoundingBox.STATUS_GRID_OK, currentVolume);
        if(!currentVolume.analyzeAngle(progressPanel, messages))
        {
            setFileStatus( StsSeismicBoundingBox.STATUS_GEOMETRY_BAD, currentVolume);
            currentVolume.progressPanelAppendTraceHeaderDescriptions(progressPanel);
            setProgressErrorStatus("Failed to analyze xy and angle geometry for " + currentVolume.getName() + ". ");
            return false;
        }
        setFileStatus( StsSeismicBoundingBox.STATUS_GEOMETRY_OK, currentVolume);
        setProgressStatus("PostStack3d " + currentVolume.getName() + " OK", StsProgressBar.INFO);
        return true;
    }
}
