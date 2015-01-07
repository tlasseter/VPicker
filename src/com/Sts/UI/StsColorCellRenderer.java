package com.Sts.UI;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.1
 */

import javax.swing.*;
import java.awt.*;

class StsColorCellRenderer extends JLabel implements ListCellRenderer
{
    private int iconWidth;
    private int iconHeight;
    private int borderWidth;
    private ColorIcon colorIcon = new ColorIcon();

    public StsColorCellRenderer()
    {
        this(16, 16, 1);
    }

    public StsColorCellRenderer(int iconWidth, int iconHeight, int borderWidth)
    {
        this.iconWidth = iconWidth;
        this.iconHeight = iconHeight;
        this.borderWidth = borderWidth;
        setOpaque(true);
        setIcon(colorIcon);
    }

    public ColorIcon createIcon(String name, Color color)
    {
        return new ColorIcon(name, color);
    }

    public Component getListCellRendererComponent( JList list, Object value,
                            int index, boolean isSelected, boolean cellHasFocus)
    {
        colorIcon.setValues((ColorIcon)value);
        setIcon(colorIcon);
        this.setText(colorIcon.name);

        if(isSelected)
        {
            setForeground(list.getSelectionForeground());
            setBackground(list.getSelectionBackground());
        }
        else
        {
            setForeground(list.getForeground());
            setBackground(list.getBackground());
        }
        return this;
    }

    class ColorIcon implements Icon
    {
        private String name;
        private Color color;

        public ColorIcon()
        {
        }

        public ColorIcon(String name, Color color)
        {
            this.name = name;
            this.color = color;
        }

        public void setValues(ColorIcon colorIcon)
        {
            if(colorIcon == null) return;
            this.name = colorIcon.name;
            this.color = colorIcon.color;
        }

        // accessors
        public Icon getIcon() { return this; }
        public String getName() { return name; }
        public Color getColor() { return color; }
        public int getIconWidth() { return iconWidth; }
        public int getIconHeight() { return iconHeight; }

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
            g.setColor(color);
            g.fillRect(x+borderWidth, y+borderWidth, iconWidth-2*borderWidth,
                    iconHeight-2*borderWidth);

            g.setColor(gColor);
        }
    }
}
