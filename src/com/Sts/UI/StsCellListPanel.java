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
import javax.swing.border.*;
import java.awt.*;

public class StsCellListPanel extends JPanel
{
    private float rowHeight;
    protected String[] items = null;

    private JScrollPane listsScrollPane = new JScrollPane();
    private Object itemPrototype;
    private ListCellRenderer cellRenderer = null;
    protected JList list = new JList();
    private GridBagLayout gridBagLayout = new GridBagLayout();

    static public final float DEFAULT_ROW_HEIGHT = 1.1f;

	public StsCellListPanel(String[] names)
	{
        this(names, DEFAULT_ROW_HEIGHT, null, null);
    }

	public StsCellListPanel(String[] names, float rowHeight, Object itemPrototype, ListCellRenderer cellRenderer)
	{
        setItems(names);
        this.rowHeight = rowHeight;
        this.itemPrototype = itemPrototype;
        this.cellRenderer = cellRenderer;
		jbInit();
	}

	private void jbInit()
	{
        setLayout(gridBagLayout);
        setBackground(Color.white);
        listsScrollPane.getViewport().add(this, null);
        int height = 0;
        if (itemPrototype != null)
        {
            if (cellRenderer != null) list.setCellRenderer(cellRenderer);
            list.setPrototypeCellValue(itemPrototype);
            height = (int)(list.getFixedCellHeight() * rowHeight);
            list.setFixedCellHeight(height);
        }
        add(list, new GridBagConstraints(0, 1, 1, 1, 0.5, 1.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                new Insets(0, 0, 0, 0), 0, 0));

        Border border = BorderFactory.createEtchedBorder();
        setBorder(border);
	}

    public void setItems(String[] names)
    {
        items = names;
        list.setListData(names);
    }

    public String[] getItems() { return items; }
    public JList getList() { return list; }

    static public void createDialog(StsCellListPanel panel)
    {
        JDialog d = new JDialog();
        d.getContentPane().add(panel);
        d.setModal(true);
        d.pack();
        d.setVisible(true);
    }

    public static void main(String[] args)
    {
        String[] names = { "Red", "Green", "Blue", "Yellow", "Red", "Green", "Blue", "Yellow", "Red", "Green", "Blue", "Yellow" };
        StsCellListPanel listPanel = new StsCellListPanel(names);
        createDialog(listPanel);
    }
}
