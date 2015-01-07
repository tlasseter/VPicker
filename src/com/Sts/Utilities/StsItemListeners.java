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
public class StsItemListeners
{
	Vector listeners = new Vector(2);

	public StsItemListeners()
	{
	}

	public synchronized void add(ItemListener listener)
	{
		if(listeners.contains(listener)) return;
		listeners.addElement(listener);
	}

	public synchronized void remove(ItemListener listener)
	{
		if (!listeners.contains(listener)) return;
		listeners.removeElement(listener);
	}

	public void fireItemStateChanged(ItemEvent e)
	{
		int count = listeners.size();
		for(int i = 0; i < count; i++)
		{
			((ItemListener)listeners.elementAt(i)).itemStateChanged(e);
		}
	}
}
