package com.Sts.Utilities;

import java.awt.event.*;
import java.util.*;

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
public class StsActionListeners
{
	Vector listeners = new Vector(2);

	public StsActionListeners()
	{
	}

	public synchronized void add(ActionListener listener)
	{
		if(listeners.contains(listener)) return;
		listeners.addElement(listener);
	}

	public synchronized void remove(ActionListener listener)
	{
		if (!listeners.contains(listener)) return;
		listeners.removeElement(listener);
	}

	public void fireActionPerformed(ActionEvent e)
	{
		int count = listeners.size();
		for(int i = 0; i < count; i++)
		{
			((ActionListener)listeners.elementAt(i)).actionPerformed(e);
		}
	}
}
