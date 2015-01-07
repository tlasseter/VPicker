package com.Sts.UI;

import com.Sts.DBTypes.*;

import javax.swing.*;
import java.awt.*;

public class StsColorListItem implements Icon
{
    protected Object object;
    protected String name;
    protected int iconWidth;
    protected int iconHeight;
    protected int borderWidth;
    private StsColor stsColor;

    public StsColorListItem(int iconWidth, int iconHeight, int borderWidth)
    {
        this(StsColor.GRAY, "gray", iconWidth, iconHeight, borderWidth);
    }

    public StsColorListItem(StsColor color, String name)
    {
        this(color, name, 16, 16, 1);
    }

    public StsColorListItem(StsColor color, String name, int iconWidth, int iconHeight, int borderWidth)
    {
        stsColor = color;
        this.name = name;
        this.iconWidth = iconWidth;
        this.iconHeight = iconHeight;
        this.borderWidth = borderWidth;
    }

    public StsColorListItem(Object object, StsColor color, String name, int iconWidth, int iconHeight, int borderWidth)
    {
        this.object = object;
        stsColor = color;
        this.name = name;
        this.iconWidth = iconWidth;
        this.iconHeight = iconHeight;
        this.borderWidth = borderWidth;
    }

    public StsColorListItem(StsMainObject object)
    {
        this(object.getStsColor(), object.getName(), 16, 16, 1);
    }

    static public StsColorListItem nullColorListItem()
    {
        return new StsColorListItem(null, StsColor.WHITE, "None", 16, 16, 1);
    }
    // accessors
    public StsColor getStsColor()
    {
        return stsColor;
    }

    public void setStsColor(StsColor c)
    {
        stsColor = c;
    }


    // accessors
    public Object getObject() { return object; }
    public void setObject(Object object) { this.object = object; }
    public Icon getIcon() { return this; }
    public void setName(String name) { this.name = name; }
    public String getName() { return name; }

    // Icon interface
    public void paintIcon(Component c, Graphics g, int x, int y)
    {
        Color gColor = g.getColor();

        // draw border
        if (borderWidth>0)
        {
            g.setColor(Color.black);  // shadow
            for (int i=0; i<borderWidth; i++)
            {
                g.drawRect(x+i, y+i, iconWidth-2*i-1, iconHeight-2*i-1);
            }
        }

        // fill icon
		if(stsColor != null)
		{
			g.setColor(stsColor.getColor());
			g.fillRect(x + borderWidth, y + borderWidth, iconWidth - 2 * borderWidth, iconHeight - 2 * borderWidth);
		}
		else
		{
			g.setColor(Color.WHITE);
			g.fillRect(x + borderWidth, y + borderWidth, iconWidth - 2 * borderWidth, iconHeight - 2 * borderWidth);
			g.setColor(Color.BLACK);
			String nullLabel = "none";
			FontMetrics fm = g.getFontMetrics();
			int stringWidth = fm.stringWidth(nullLabel);
			int stringHeight = fm.getHeight();
			g.drawString("none", x + (iconWidth - stringWidth)/2, y + (iconHeight - stringHeight)/2);
		}
		g.setColor(gColor);
    }
    public int getIconWidth() { return iconWidth; }
    public int getIconHeight() { return iconHeight; }

	static public StsColor[] getColorsFromItems(StsColorListItem[] items)
	{
		if(items == null) return new StsColor[0];
		int nColors = items.length;
		StsColor[] colors = new StsColor[nColors];
		for(int n = 0; n < nColors; n++)
			colors[n] = items[n].stsColor;
		return colors;
	}

    public String toString()
	{
		return name;
	}

    public boolean equals(Object object)
    {
        if(object == null) return false;
		if(name != null)
        	return toString().equals(object.toString());
		else
			return stsColor.equals(((StsColorListItem)object).stsColor);
    }
}
