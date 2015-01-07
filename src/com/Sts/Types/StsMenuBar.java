
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import java.awt.*;

public abstract class StsMenuBar extends MenuBar
{
//	public StsMenuBar() { }

	public abstract void setEnabledMenuItems(boolean enable);
	public abstract void setCheckboxUseDisplayLists(boolean enable);
	public abstract void setCheckboxView3dCursors(boolean enable);
}
