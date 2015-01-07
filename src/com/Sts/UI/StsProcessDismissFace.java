package com.Sts.UI;

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
public interface StsProcessDismissFace
{
	static public final int PROCESS = 1;
	static public final int DISMISS = 2;
	static public final int CLOSING = 3;

	public void process();
	public void dismiss();
}
