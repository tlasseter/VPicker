package com.Sts.UI.Beans;

import com.Sts.DBTypes.StsColor;
import com.Sts.DBTypes.StsSpectrum;
import com.Sts.UI.StsButton;
import com.Sts.UI.StsColorListItem;
import com.Sts.UI.StsColorListRenderer;
import com.Sts.Utilities.StsException;
import com.Sts.Utilities.StsParameters;
import com.Sts.Utilities.StsToolkit;

import java.awt.*;

public class StsColorItemBean extends StsJPanel
{
	public StsColorItemBean()
    {
    }

	public StsColorItemBean(StsStringFieldBean nameBean, StsComboBoxFieldBean colorBean)
    {
		addToRow(nameBean);
		addEndRow(colorBean);
    }
}

