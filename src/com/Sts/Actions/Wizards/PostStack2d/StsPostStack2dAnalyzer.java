package com.Sts.Actions.Wizards.PostStack2d;

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

public class StsPostStack2dAnalyzer extends StsPostStackAnalyzer
{
    static StsPostStack2dAnalyzer analyzer = null;

    static public StsPostStack2dAnalyzer getAnalyzer(StsPostStack2dWizard wizard, StsProgressPanel progressPanel, StsTablePanelNew volumeStatusPanel)
    {
        if(analyzer != null)
         {
             // cancel any running analyzer; analyze routines are synchronized,
             // so new analyzer will wait until any old one is finished
             analyzer.cancelProcess();
         }
         else
             analyzer = new StsPostStack2dAnalyzer();

        analyzer.initialize(wizard, progressPanel, volumeStatusPanel);
        return analyzer;
    }

    public boolean analyzeGrid(StsProgressPanel progressPanel, boolean messages, StsSeismicBoundingBox volume)
    {
        StsSegyLine2d currentVolume = (StsSegyLine2d)volume;
        if(!currentVolume.analyzeGrid(progressPanel, messages))
        {
            setFileStatus( StsSeismicBoundingBox.STATUS_GEOMETRY_BAD, currentVolume);
            setProgressErrorStatus("Failed to analyze 2d line geometry for " + currentVolume.getName() + ". ");
            return false;
        }
        setFileStatus( StsSeismicBoundingBox.STATUS_GEOMETRY_OK, currentVolume);
        setProgressStatus("PostStack2d " + currentVolume.getName() + " OK", StsProgressBar.INFO);
        return true;
    }

    public boolean analyzeAngle(StsProgressPanel progressPanel, boolean messages, StsSeismicBoundingBox volume)
    {
        return true;
    }
}