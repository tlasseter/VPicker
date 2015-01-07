
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import java.awt.*;

public class DrawnRectangle extends Rectangle
{
    protected static int _defaultThickness = 2;
    private int thick;
    private Color lineColor, fillColor;

    public DrawnRectangle()
    {
        this(_defaultThickness, 0, 0, 0, 0);
    }
    public DrawnRectangle(int thick)
    {
        this(thick, 0, 0, 0, 0);
    }
    public DrawnRectangle(int x, int y, int w, int h)
    {
        this(_defaultThickness, x, y, w, h);
    }
    public DrawnRectangle(int thick, int x, int y, int w, int h)
    {
        this.thick    = thick;
        setBounds(x,y,w,h);
    }
    public int  getThickness() { return thick;}
    public void setThickness(int thick) { this.thick = thick; }

    public void setLineColor(Color lineColor)
    {
        this.lineColor = lineColor;
    }
    public void setFillColor(Color fillColor)
    {
        this.fillColor = fillColor;
    }
    public void fill(Graphics g)
    {
        fill(g, getFillColor());
    }
    public Color getLineColor()
    {
        if(lineColor == null)
            lineColor = SystemColor.controlShadow;
        return lineColor;
    }
    public Color getFillColor()
    {
        if(fillColor == null)
            fillColor = Color.cyan;
        return fillColor;
    }
    public Rectangle getInnerBounds()
    {
        return new Rectangle(x+thick, y+thick,
                             width-(thick*2), height-(thick*2));
    }
    public void paint(Graphics g)
    {
        paintFlat(g, getLineColor());
    }

    private void paintFlat(Graphics g, Color color)
    {
        if(g != null)
        {
            g.setColor(color);
            for(int i=0; i < thick; ++i)
                g.drawRect(x+i, y+i, width-(i*2)-1, height-(i*2)-1);
        }
    }
    public void fill(Graphics g, Color color)
    {
        if(g != null)
        {
            Rectangle r = getInnerBounds();
            g.setColor(color);
            g.fillRect(r.x, r.y, r.width, r.height);
            setFillColor(color);
        }
    }
    protected Color brighter()
    {
        return getLineColor().brighter().brighter().brighter().brighter();
    }
}
