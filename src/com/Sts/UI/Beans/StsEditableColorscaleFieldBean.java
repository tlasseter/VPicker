package com.Sts.UI.Beans;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.Actions.Wizards.Color.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;
import java.awt.event.*;

/** A fieldBean used currently in the objectPanel to display an editable colorscale.
 *  The bean objectType is an StsColorscale, i.e., we get/set colorscales into this
 *  bean as the panel holding the bean changes its panelObject as when we select different
 *  instances of an class which displays this fieldbean.
 */

public class StsEditableColorscaleFieldBean extends StsFieldBean
{
	StsColorscalePanel colorPanel = new StsColorscalePanel(true, StsColorscalePanel.COLORSCALE); // Setting mode to Colorscale (0)
	private StsButton editButton = new StsButton("  Edit Colorscale  ", "Edit colorscale colors, range, and style.", this, "doColorscaleAction");

    public StsEditableColorscaleFieldBean()
    {
		initialize();
	}

    public StsEditableColorscaleFieldBean(Class c, String fieldName)
    {
        super.classInitialize(c, fieldName, null);
		initialize();
    }

    public StsEditableColorscaleFieldBean(Object instance, String fieldName)
    {
        initialize(instance, fieldName);
		initialize();
    }

	private void initialize()
	{
//		addVerticalSpacer = true;
		add(colorPanel);
 		add(editButton);
		colorPanel.setLabelsOn(true);
		this.gridwidth = 2;
	}

    public Component[] getBeanComponents() { return new Component[] { this }; }
//    public StsEditableColorscalePanel getPanel() { return panel; }

    public Object getValueObject() { return colorPanel.getColorscale(); }

    public void doSetValueObject(Object object)
    {
        if(object != null && !(object instanceof StsColorscale))
        {
            super.outputSetValueObjectException(object, "StsColorscale");
        }
        colorPanel.setColorscale((StsColorscale)object);
        colorPanel.repaint();
    }

	public String toString() {return StsParameters.NONE_STRING; }
	public Object fromString(String string) { return null; }

    public void actionPerformed(ActionEvent e)
    {
        setValueInPanelObjects();
    }

	public void setHistogram(float[] dataHist_)
	{
        if(colorPanel == null) return;
        final float[] dataHist = dataHist_;
        StsToolkit.runLaterOnEventThread
        (
            new Runnable()
            {
                public void run()
                {
                    colorPanel.setHistogram(dataHist);
                }
            }
        );
	}
    
    public float[] getHistogram() { return colorPanel.getHistogram(); }

	public void doColorscaleAction()
	{
		StsColorscale colorscale = colorPanel.getColorscale();
//		float[] data = getHistogram();
		StsModel model = StsObject.getCurrentModel();
		if(model == null) // this is for standalone testing
		{
			StsActionManager actionManager = new StsActionManager(model);
			StsColorscaleAction action = new StsColorscaleAction(actionManager, colorPanel);
			action.start();
		}
		else
		{
			StsActionManager actionManager = model.win3d.getActionManager();
			actionManager.startAction(StsColorscaleAction.class, new Object[] { colorPanel} );
		}
	}

	public void setColorscale(StsColorscale colorscale)
	{
		colorPanel.setLabelsOn(true);
		editButton.setVisible(colorscale != null);
		colorPanel.setColorscale(colorscale);
	}

	public StsColorscale getColorscale()
	{
		return colorPanel.getColorscale();
	}

	public void setRange(float[] range)
	{
		colorPanel.setRange(range);
	}

	public void setRange(float newMin, float newMax)
	{
		colorPanel.setRange(newMin, newMax);
	}

    public StsEditableColorscaleFieldBean copy(Object beanObject)
	{
		try
		{
			StsEditableColorscaleFieldBean beanCopy = new StsEditableColorscaleFieldBean();
            beanCopy.initialize(beanObject, fieldName, editable, getLabelString());
			return beanCopy;
		}
		catch(Exception e)
		{
			return null;
		}
	}

    // This main() doesn't allow for editing because an StsColorscaleAction is not being used.  TJL 3/18/07
    static public void main(String[] args)
	{
		try
		{
			StsModel model = new StsModel();
			StsActionManager actionManager = new StsActionManager(model);
			CFBTest test = new CFBTest();
			StsEditableColorscaleFieldBean colorscaleBean = new StsEditableColorscaleFieldBean();
			colorscaleBean.initialize(test, "colorscale");
			StsObjectPanel panel = StsObjectPanel.constructEmptyPanel();
			panel.add(colorscaleBean);
			StsToolkit.createDialog(panel, true);
		}
		catch(Exception e)
		{
			e.printStackTrace();
        }
	}
}
