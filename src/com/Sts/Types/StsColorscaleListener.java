
//Title:        S2S: Seismic-to-simulation: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System


package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;

import javax.swing.event.*;
import java.awt.event.*;

public class StsColorscaleListener implements ActionListener, ChangeListener, ItemListener
{
    StsColorscalePanel colorscalePanel = null;

    public StsColorscaleListener(StsColorscalePanel p)
    {
        this.colorscalePanel = p;
//        p.getColorscale().addItemListener(this);
    }

    public void actionPerformed(ActionEvent e)
    {
        StsMessageFiles.infoMessage("actionPerformed Performed...");
        return;
    }

    public void stateChanged(ChangeEvent e)
    {
        StsMessageFiles.infoMessage("stateChanged Performed...");
        return;
    }

    public void itemStateChanged(ItemEvent e)
    {
        StsMessageFiles.infoMessage("itemStateChanged Performed...");
        if( e.getItem() instanceof StsColorscale )
        {
            StsColorscale csp = (StsColorscale) e.getItem();
            updateSeismicColors(csp);
        }
        return;
    }

    private void updateSeismicColors(StsColorscale cs)
    {
        StsModel model = StsSerialize.getCurrentModel();
        StsSeismicVolume seismicVolume = (StsSeismicVolume)model.getStsClass(StsSeismicVolume.class).getCurrentObject();
        if(seismicVolume == null) return;
        seismicVolume.resetColors();
        model.win3dDisplayAll();
    }
}
