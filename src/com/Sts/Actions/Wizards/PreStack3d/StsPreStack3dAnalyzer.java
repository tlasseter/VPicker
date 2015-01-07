package com.Sts.Actions.Wizards.PreStack3d;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Types.PreStack.*;
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

public class StsPreStack3dAnalyzer extends StsPostStackAnalyzer
{
    private static StsPreStack3dAnalyzer analyzer;

    static public StsPreStack3dAnalyzer getAnalyzer(StsSeismicWizard wizard, StsProgressPanel progressPanel, StsTablePanelNew volumeStatusPanel)
    {
        if(analyzer != null)
        {
            // cancel any running analyzer; analyze routines are synchronized,
            // so new analyzer will wait until any old one is finished
            analyzer.clearProcess();
        }
        else
            analyzer = new StsPreStack3dAnalyzer();
        analyzer.initialize(wizard, progressPanel, volumeStatusPanel);
        return analyzer;
    }

    public boolean analyzeGrid(StsProgressPanel progressPanel, boolean messages, StsSeismicBoundingBox volume)
    {
        StsPreStackSegyLine currentVolume = (StsPreStackSegyLine)volume;
        if(!currentVolume.analyzeGeometry())
        {
            setFileStatus( StsSeismicBoundingBox.STATUS_GEOMETRY_BAD, currentVolume);
//            currentVolume.progressPanelAppendTraceHeaderDescriptions(progressPanel);
            setProgressErrorStatus("Failed to analyze xy and angle geometry for " + currentVolume.getName() + ". ");
            return false;
        }
        setFileStatus( StsSeismicBoundingBox.STATUS_GEOMETRY_OK, currentVolume);
        setProgressStatus("PostStack3d " + currentVolume.getName() + " OK", StsProgressBar.INFO);
        return true;
    }
}
