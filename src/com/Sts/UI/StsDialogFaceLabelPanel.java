package com.Sts.UI;

import com.Sts.Interfaces.StsDialogFace;
import com.Sts.UI.Beans.StsGroupBox;
import com.Sts.UI.Beans.StsStringFieldBean;
import com.Sts.Utilities.StsToolkit;

import javax.swing.*;
import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: May 6, 2009
 * Time: 5:23:25 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsDialogFaceLabelPanel implements StsDialogFace
{
    String panelname;
    public String string = "";
	JLabel label;
    boolean okSelected = false;

    public StsDialogFaceLabelPanel(String panelname, String stringLabel)
	{
        this.panelname = panelname;
        label = new JLabel(stringLabel);
    }

	public void dialogSelectionType(int type)
	{
        okSelected = (type == StsDialogFace.OK);
    }
    public boolean okSelected() { return okSelected; }
	public Component getPanel(boolean val) { return getPanel(); }

    public Component getPanel()
	{
		StsGroupBox groupBox = new StsGroupBox(panelname);
		groupBox.add(label);
		return groupBox;
	}

    public StsDialogFace getEditableCopy()
    {
        return (StsDialogFace) StsToolkit.copyObjectNonTransientFields(this);
    }
}