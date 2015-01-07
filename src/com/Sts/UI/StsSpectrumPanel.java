package com.Sts.UI;

import com.Sts.DBTypes.*;

import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsSpectrumPanel extends JPanel
{
    private JPanel panel1 = new JPanel();
    StsColorscalePanel colorPanel = new StsColorscalePanel(true, StsColorscalePanel.SPECTRUM); // Setting mode to Spectrum (1)
    private BorderLayout borderLayout1 = new BorderLayout();
    private Border border1;

    public StsSpectrumPanel()
    {
        try
        {
            jbInit();
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

    private void jbInit() throws Exception
    {
        border1 = BorderFactory.createCompoundBorder(BorderFactory.createBevelBorder(BevelBorder.LOWERED,Color.white,Color.white,new Color(148, 145, 140),new Color(103, 101, 98)),BorderFactory.createEmptyBorder(4,4,4,4));
//       this.setLayout(verticalLayout);

        this.setMinimumSize(new Dimension(55, 350));
        this.setPreferredSize(new Dimension(55, 350));
        this.add(colorPanel, null);
        panel1.setLayout(borderLayout1);
        this.add(panel1, null);
    }

    public void setSpectrum(StsSpectrum spectrum)
    {
        colorPanel.setLabelsOn(true);
        colorPanel.setSpectrum(spectrum);
    }
}
