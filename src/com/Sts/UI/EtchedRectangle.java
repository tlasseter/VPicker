
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import java.awt.*;

public class EtchedRectangle extends DrawnRectangle
{
    public static final int IN = 1;
	public static final int OUT = 1;
    protected static int _defaultEtching = IN;
    private int etching;

    public EtchedRectangle()
    {
        this(_defaultEtching, _defaultThickness, 0, 0, 0, 0);
    }
    public EtchedRectangle(int thickness)
    {
        this(_defaultEtching, thickness, 0, 0, 0, 0);
    }
    public EtchedRectangle(int x, int y, int w, int h)
    {
        this(_defaultEtching, _defaultThickness, x, y, w, h);
    }
    public EtchedRectangle(int thickness, int x, int y, int w, int h)
    {
        this(_defaultEtching, thickness, x, y, w, h);
    }
    public EtchedRectangle(int etching, int thickness, int x, int y, int w, int h)
    {
        super(thickness, x, y, w, h);
        this.etching = etching;
    }
    public void    etchedIn  () { etching = IN;        }
    public void    etchedOut () { etching = OUT;       }
    public boolean isEtchedIn() { return etching == IN;}

    public void paint(Graphics g)
    {
        if(etching == IN) paintEtchedIn(g);
        else paintEtchedOut(g);
    }
    public void paintEtchedIn(Graphics g)
    {
        if(g != null)
            paintEtched(g, getLineColor(), brighter());

        etchedIn();
    }
    public void paintEtchedOut(Graphics g)
    {
        if(g != null)
            paintEtched(g, brighter(), getLineColor());

        etchedOut();
    }
    private void paintEtched(Graphics g, Color topLeft, Color bottomRight)
    {
        int  thickness = getThickness();
        int  w = width  - thickness;
        int  h = height - thickness;

        g.setColor(topLeft);
        for(int i=0; i < thickness/2; ++i)
            g.drawRect(x+i, y+i, w, h);

        g.setColor(bottomRight);

        for(int i=0; i < thickness/2; ++i)
            g.drawRect(x+(thickness/2)+i,
                       y+(thickness/2)+i, w, h);

    }
}
