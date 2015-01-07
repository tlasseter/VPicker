package com.Sts.UI.Histogram;

import com.Sts.Types.*;
import com.Sts.UI.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Mar 26, 2008
 * Time: 11:27:57 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsSeismicVolumesHistogramPanel extends StsHistogramPanel
{

    StsSeismicBoundingBox[] selectedVolumes;
    private StsTablePanelNew volumeStatusPanel = null;

    public StsSeismicVolumesHistogramPanel(StsTablePanelNew volumeStatusPanel)
    {
        super(HORIZONTAL, true);
        this.volumeStatusPanel = volumeStatusPanel;
    }

    public void setSelectedVolumes(StsSeismicBoundingBox[] selectedVolumes)
    {
        this.selectedVolumes = selectedVolumes;
    }

    public void setClipDataMin(float clipDataMin)
    {
        super.setClipDataMin(clipDataMin);
        if( selectedVolumes == null) return;
        for(int n = 0; n < selectedVolumes.length; n++)
            selectedVolumes[n].setDataMin(clipDataMin);
        volumeStatusPanel.repaint();
    }

    public void setClipDataMax(float clipDataMax)
    {
        super.setClipDataMax(clipDataMax);
        if( selectedVolumes == null) return;
        for(int n = 0; n < selectedVolumes.length; n++)
            selectedVolumes[n].setDataMax(clipDataMax);
        volumeStatusPanel.repaint();
    }
}
