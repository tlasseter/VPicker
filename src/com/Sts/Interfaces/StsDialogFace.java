package com.Sts.Interfaces;

import java.awt.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */
public interface StsDialogFace
{
	static public final int OK = 1;
	static public final int APPLY = 2;
	static public final int CANCEL = 3;
	static public final int CLOSING = 4;
	static public final int PROCESS = 5;
	static public final int DISMISS = 6;

	public void dialogSelectionType(int type);
	public Component getPanel();
	public Component getPanel(boolean expand);
    public StsDialogFace getEditableCopy();
}
