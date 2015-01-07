package com.Sts.UI;

import com.Sts.Actions.Wizards.Color.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;

import javax.swing.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */


public class StsEditableColorscalePanel extends StsJPanel
{
    private StsButton editButton = new StsButton(" Edit colorscale ", "Change colors, ranges, and styles.", this, "doEditColorscale");
    private JPanel panel1 = new JPanel();
    StsColorscalePanel colorPanel = new StsColorscalePanel(true, StsColorscalePanel.COLORSCALE); // Setting mode to Colorscale (0)

    public StsEditableColorscalePanel()
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
		add(colorPanel);
		add(editButton);
    }

	public void doEditColorscale()
	{
		StsColorscale colorscale = colorPanel.getColorscale();
		float[] data = colorPanel.getHistogram();
		StsModel model = StsSerialize.getCurrentModel();
		if(data != null)
		{
			model.mainWindowActionManager.startAction(StsColorscaleAction.class, new Object[]
			{
				colorscale, data
			});
		}
		else
		{
			model.mainWindowActionManager.startAction(StsColorscaleAction.class, new Object[]
			{
				colorscale
			});
		}
	}

    public void setColorscale(StsColorscale colorscale)
    {
        colorPanel.setLabelsOn(true);
        editButton.setVisible(colorscale != null);
        colorPanel.setColorscale(colorscale);
    }

    public StsColorscale getColorscale() { return colorPanel.getColorscale(); }

    public void setRange(float[] range)
    {
        colorPanel.setRange(range);
    }
    public void setRange(float newMin, float newMax)
    {
        colorPanel.setRange(newMin, newMax);
    }
}
