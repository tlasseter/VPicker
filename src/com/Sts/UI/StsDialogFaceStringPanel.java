package com.Sts.UI;

import com.Sts.Interfaces.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: May 6, 2009
 * Time: 5:23:25 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsDialogFaceStringPanel implements StsDialogFace
{
    String panelname;
    String initialString;
    public String string = "";
	StsStringFieldBean stringBean;

    public StsDialogFaceStringPanel(String panelname, String stringLabel, String initialString)
	{
        string = initialString;
        this.initialString = initialString;
        this.panelname = panelname;
        stringBean = new StsStringFieldBean(this, "string", stringLabel);
    }

	public void dialogSelectionType(int type)
	{
        boolean okSelected = (type == StsDialogFace.OK);
        if(!okSelected) string = initialString;

    }
	public Component getPanel(boolean val) { return getPanel(); }

    public Component getPanel()
	{
		StsGroupBox groupBox = new StsGroupBox(panelname);
		groupBox.add(stringBean);
		return groupBox;
	}

    public StsDialogFace getEditableCopy()
    {
        return (StsDialogFace) StsToolkit.copyObjectNonTransientFields(this);
    }

    public void setString(String s) { string = s; }
	public String getString() { return string; }
}
