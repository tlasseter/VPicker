
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

// dialog box to set a text field

package com.Sts.UI;

import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;

public class StsColorListComboBox extends JComboBox
{
    private int iconHeight;
    private int iconWidth;
    private ColorCellRenderer renderer;

	public StsColorListComboBox()
	{
        this(16, 16, false);
    }

	public StsColorListComboBox(int iconWidth, int iconHeight, boolean hasText)
	{
        this((String)null, null, iconWidth, iconHeight, hasText);
    }

	public StsColorListComboBox(String name, StsColor color, int iconWidth, int iconHeight, boolean hasText)
	{
        super();
//        setLightWeightPopupEnabled(false);
        this.setBackground(Color.white);
        this.setForeground(Color.black);
        renderer = new ColorCellRenderer(iconWidth, iconHeight, hasText);
        setRenderer(renderer);
        addItem(name, color);
	}

	public StsColorListComboBox(String[] names, StsColor[] colors, int iconWidth, int iconHeight, boolean hasText)
	{
        super();
//        setLightWeightPopupEnabled(false);
        this.setBackground(Color.white);
        this.setForeground(Color.black);
        renderer = new ColorCellRenderer(iconWidth, iconHeight, hasText);
        setRenderer(renderer);
        setItems(names, colors);
	}

    public void addItem(String name, StsColor color)
    {
        ColorCellRenderer.Item item = renderer.createItem(name, color);
        this.addItem(item);
    }

    public StsColor getSelectedColor()
    {
        ColorCellRenderer.Item item = (ColorCellRenderer.Item)getSelectedItem();
        return item.getColor();
    }

    public String getSelectedName()
    {
        ColorCellRenderer.Item item = (ColorCellRenderer.Item)getSelectedItem();
        return item.getName();
    }

    public void setItems(String[] names, StsColor[] colors)
    {
        int nItems;

        if (colors==null && names==null) return;
        if(names == null)
        {
            nItems = colors.length;
            for (int i=0; i<nItems; i++)
            {
                ColorCellRenderer.Item item = renderer.createItem("", colors[i]);
                addItem(item);
            }
        }
        else if(colors == null)
        {
            nItems = names.length;
            for (int i=0; i<nItems; i++)
            {
                ColorCellRenderer.Item item = renderer.createItem(names[i], null);
                addItem(item);
            }
        }
        else
        {
            nItems = Math.min(colors.length, names.length);
            for (int i=0; i<nItems; i++)
            {
                ColorCellRenderer.Item item = renderer.createItem(names[i], colors[i]);
                addItem(item);
            }
        }
    }

    static public StsColorListComboBox createTest()
    {
        String[] names = StsColor.colorNames8;
        StsColor[] colors = StsColor.colors8;
        StsColorListComboBox comboBox = new StsColorListComboBox(names, colors, 30, 15, true);
//        comboBox.setLightWeightPopupEnabled(false);
        return comboBox;
    }

    public static void main(String[] args)
    {
        StsColorListComboBox comboBox = StsColorListComboBox.createTest();
        StsToolkit.createDialog(comboBox);
    }
}

class ColorCellRenderer extends JPanel implements ListCellRenderer
{
    private int iconWidth;
    private int iconHeight;
    private boolean hasText;

    private Item item = new Item();

    JLabel colorLabel = new JLabel();
    JLabel textLabel = new JLabel();
    GridBagLayout gridBagLayout = new GridBagLayout();

    private Border
        blackBorder = BorderFactory.createLineBorder(Color.black, 2),
        emptyBorder = BorderFactory.createEmptyBorder(2, 2, 2, 2);

    public ColorCellRenderer(boolean hasText)
    {
        this(16, 16, hasText);
    }

    public ColorCellRenderer(int iconWidth, int iconHeight, boolean hasText)
    {
        this.hasText = hasText;

//        this.setLayout(gridBagLayout);
        this.setLayout(new BorderLayout());
        this.setOpaque(true);
        colorLabel.setOpaque(true);
        colorLabel.setPreferredSize(new Dimension(iconWidth, iconHeight));
        colorLabel.setBorder(BorderFactory.createLineBorder(Color.black));
        this.add(colorLabel, BorderLayout.WEST);
    /*
        this.add(colorLabel,      new GridBagConstraints(0, 0, 1, 1, 0.5, 0.5
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(4, 4, 4, 4), 0, 0));
    */
        if(hasText)
        {
            textLabel.setFont(new java.awt.Font("Dialog", 0, 12));
//            textLabel.setForeground(Color.black);
            textLabel.setPreferredSize(new Dimension(80, 15));
            textLabel.setHorizontalAlignment(SwingConstants.LEFT);
            this.add(textLabel, BorderLayout.EAST);
        /*
            this.add(textLabel,      new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0
                ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(4, 4, 4, 4), 0, 0));
        */
        }
    }

    public Component getListCellRendererComponent( JList list, Object value,
                            int index, boolean isSelected, boolean cellHasFocus)
    {
        item.setValues((Item)value);
        colorLabel.setBackground(item.color.getColor());
        if(hasText) textLabel.setText(item.name);

        if(isSelected)
        {
            this.setBorder(blackBorder);
            setBackground(Color.white);
            textLabel.setBackground(Color.white);
//            setForeground(list.getSelectionForeground());
//            setBackground(list.getSelectionBackground());
        }
        else
        {
            this.setBorder(emptyBorder);
            setBackground(Color.white);
//            setForeground(list.getForeground());
//            setBackground(list.getBackground());
        }
        return this;
    }

    public Item createItem(String name, StsColor color)
    {
        return new Item(name, color);
    }

    class Item
    {
        private String name;
        private StsColor color;

        public Item()
        {
            this("gray", StsColor.GRAY);
        }

        public Item(String name, StsColor color)
        {
            this.name = name;
            this.color = color;
        }

        public void setValues(Item item)
        {
            if(item == null) return;
            this.name = item.name;
            this.color = item.color;
        }

        // accessors
        public String getName() { return name; }
        public StsColor getColor() { return color; }
    }
}
