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

public class StsCellListComboBox extends JComboBox
{
    protected String[] items = null;

    private ListCellRenderer cellRenderer;

	public StsCellListComboBox(String[] names)
	{
        this(names, null);
        this.setKeySelectionManager(new NullKeySelectionManager()); // turns off inadvertent key selection
    }

	public StsCellListComboBox(String[] names, ListCellRenderer cellRenderer)
	{
        setItems(names);
        setLightWeightPopupEnabled(false);
        this.cellRenderer = cellRenderer;
        this.setSelectedIndex(0);
		jbInit();
	}

	private void jbInit()
	{
        if (cellRenderer != null) this.setRenderer(cellRenderer);
	}

    public void setItems(String[] names)
    {
        items = names;
        for(int n = 0; n < names.length; n++)
            super.addItem(names[n]);
    }

    public String[] getItems() { return items; }

    static public void createDialog(StsCellListComboBox comboBox)
    {
        JDialog d = new JDialog();
        d.getContentPane().add(comboBox);
        d.setModal(true);
        d.pack();
        d.setVisible(true);
    }

    static public StsCellListComboBox createTest()
    {
        String[] names = { "Red", "Green", "Blue", "Yellow", "Red", "Green", "Blue", "Yellow", "Red", "Green", "Blue", "Yellow" };
        StsCellListComboBox comboBox = new StsCellListComboBox(names);
        comboBox.setLightWeightPopupEnabled(false);
        return comboBox;
    }

    public static void main(String[] args)
    {
        StsCellListComboBox comboBox = StsCellListComboBox.createTest();
        createDialog(comboBox);
    }
}

class NullKeySelectionManager implements JComboBox.KeySelectionManager
{
    public NullKeySelectionManager()
    {
    }

    public int selectionForKey(char aKey, ComboBoxModel aModel)
    {
        return -1;
    }
}
