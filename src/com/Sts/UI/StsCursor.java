
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.image.*;

public class StsCursor
{
    private Component component;
    private Cursor oldCursor = Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR);
	private int cursorType = Cursor.DEFAULT_CURSOR;
	static public final int DEFAULT = Cursor.DEFAULT_CURSOR; // 0
	// these enumerations need to match mouse modes
    static public final int ZOOM = 1;
    static public final int PAN = 2;
    static public final int ROTATE = 3;
	static public final int RECTZOOM = 4;
	static public final int PANVERTICAL= 5;

	static public final String[] mouseModeStrings = new String[] { "DEFAULT", "ZOOM", "PAN", "ROTATE", "RECTZOOM", "PANVERTICAL"};
	static Cursor[] stsCursors;

    /** constructor with no change of cursor */
    public StsCursor(Component component)
    {
        if (component == null)
            StsException.outputException(new StsException(StsException.WARNING,
                "StsCursor.StsCursor:  unable to build cursor for null component."));
        else
        {
            this.component = component;
            oldCursor = component.getCursor();
        }
		if(stsCursors == null) constructStsCursors();
    }

    /** constructor with cursor change */
    public StsCursor(Component component, int cursorType)
    {
        if (component==null)
            StsException.systemError("StsCursor.StsCursor: unable to build cursor for null component.");
        else
        {
            this.component = component;
			this.cursorType = cursorType;
            setCursor(cursorType);
        }
		if(stsCursors == null) constructStsCursors();
    }

    /** set the cursor for the component */
    public void setCursor(int cursorType)
    {
		if(this.cursorType == cursorType) return;
        if (!cursorTypeIsValid(cursorType))
            StsException.outputException(new StsException(StsException.WARNING,
                "StsCursor.setCursor:  unknown cursor type requested."));
        else
        {
            oldCursor = component.getCursor();
            component.setCursor(Cursor.getPredefinedCursor(cursorType));
        }
    }

    public void setStsCursor(int cursorType)
    {
		if(this.cursorType == cursorType) return;
        oldCursor = component.getCursor();
        component.setCursor(stsCursors[cursorType]);
    }

	public int getStsCursor()
	{
		return this.cursorType;
    }

    /** restore the previous cursor */
    public void restoreCursor()
    {
		cursorType = Cursor.DEFAULT_CURSOR;
        if(component == null) return;
        oldCursor = Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR);
		component.setCursor(oldCursor);
   }

    public Cursor getCursor()
    {
        return component.getCursor();
    }

    /** check for valid cursor */
    static public boolean cursorTypeIsValid(int cursorType)
    {
        switch (cursorType)
        {
            case Cursor.CROSSHAIR_CURSOR:
            case Cursor.DEFAULT_CURSOR:
            case Cursor.E_RESIZE_CURSOR:
            case Cursor.HAND_CURSOR:
            case Cursor.MOVE_CURSOR:
            case Cursor.N_RESIZE_CURSOR:
            case Cursor.NE_RESIZE_CURSOR:
            case Cursor.NW_RESIZE_CURSOR:
            case Cursor.S_RESIZE_CURSOR:
            case Cursor.SE_RESIZE_CURSOR:
            case Cursor.SW_RESIZE_CURSOR:
            case Cursor.TEXT_CURSOR:
            case Cursor.W_RESIZE_CURSOR:
            case Cursor.WAIT_CURSOR:
                return true;
        }
        return false;
    }


	static private void constructStsCursors()
	{
		java.net.URL url;
		Toolkit tk = Toolkit.getDefaultToolkit();
		Image img;

		stsCursors = new Cursor[6];

		try
		{
			stsCursors[DEFAULT] = Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR);

			url = StsCursor.class.getResource("Icons/zoomCursor.gif");
			img = tk.createImage( (ImageProducer)url.getContent());
			stsCursors[ZOOM] = tk.createCustomCursor(img, new Point(16, 16), "zoom");

	        //url = StsCursor.class.getResource("Icons/zoomCursor.gif");
			//img = tk.createImage( (ImageProducer)url.getContent());
			stsCursors[RECTZOOM] = Cursor.getPredefinedCursor(Cursor.NW_RESIZE_CURSOR);


			url = StsCursor.class.getResource("Icons/panCursor.gif");
			img = tk.createImage( (ImageProducer)url.getContent());
			stsCursors[PAN] = tk.createCustomCursor(img, new Point(16, 16), "pan");

			url = StsCursor.class.getResource("Icons/rotateCursor.gif");
			img = tk.createImage( (ImageProducer)url.getContent());
			stsCursors[ROTATE] = tk.createCustomCursor(img, new Point(16, 16), "rotate");

	         stsCursors[PANVERTICAL] = new Cursor(Cursor.HAND_CURSOR);
		}
		catch (Exception e) { e.printStackTrace(); }
    }

    static void main(String[] args)
    {
        try
        {
            JPanel panel = new JPanel();
            panel.setSize(new Dimension(200,200));
    //        Toolkit tk = panel.getToolkit();

            java.net.URL url = StsCursor.class.getResource("zoomCursor.gif");
            Toolkit tk = Toolkit.getDefaultToolkit();
            Image img = tk.createImage((ImageProducer)url.getContent());
            Cursor cursor = tk.createCustomCursor(img, new Point(16, 16), "zoom");
            panel.setCursor(cursor);

            JFrame frame = new JFrame();
            frame.getContentPane().add(panel);
            frame.pack();
            frame.setVisible(true);
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }
}
