

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import java.awt.event.*;

public class StsMouseButton
{
	public StsMousePoint mousePoint = new StsMousePoint();
    public StsMousePoint mouseDelta = new StsMousePoint();
    public int  actionID;

	public StsMouseButton()
	{
	}

    public void setState(MouseEvent e)
    {
    	actionID = e.getID();

     	if(actionID == MouseEvent.MOUSE_DRAGGED)
        {
    		int newX = e.getX();
            int newY = e.getY();

            mouseDelta.x = newX - mousePoint.x;
            mouseDelta.y = -newY + mousePoint.y;

            mousePoint.x = newX;
            mousePoint.y = newY;
        }
        else
        {
        	mouseDelta.x = 0;
            mouseDelta.y = 0;
        	mousePoint.x = e.getX();
        	mousePoint.y = e.getY();
        }
    }

    public void setStateAndID(MouseEvent e, int ID)
    {
    	actionID = ID;

     	if(actionID == MouseEvent.MOUSE_DRAGGED)
        {
    		int newX = e.getX();
            int newY = e.getY();

            mouseDelta.x = newX - mousePoint.x;
            mouseDelta.y = -newY + mousePoint.y;

            mousePoint.x = newX;
            mousePoint.y = newY;
        }
        else
        {
        	mouseDelta.x = 0;
            mouseDelta.y = 0;
        	mousePoint.x = e.getX();
        	mousePoint.y = e.getY();
        }
    }

    public StsMousePoint getMousePoint()
    {
    	return mousePoint;
    }

    public StsMousePoint getMouseDelta()
    {
    	return mouseDelta;
    }
}
