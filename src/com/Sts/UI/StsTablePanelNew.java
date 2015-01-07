package com.Sts.UI;

import com.Sts.DBTypes.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Table.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import java.awt.*;
import java.util.*;

public class StsTablePanelNew extends StsGroupBox implements ListSelectionListener, TableModelListener
{
    private StsEditableTableModel tableModel;
    private JTable table;
    private int scrollPanelWidth = defaultWidth;
    private int scrollPaneHeight = defaultHeight;
    private transient Vector listSelectionListeners = null;
    private transient Vector tableModelListeners = null;
    private java.util.List selectRowNotifyListeners = new LinkedList();
    private java.util.List deleteRowNotifyListeners = new LinkedList();
    private Color selectionColor = null;
    public static final byte NOT_HIGHLIGHTED = 0; // White
    public static final byte HIGHLIGHTED = 1; // Light Blue
    public static final byte NOT_EDITABLE = 2; // Grey
    public static final byte SELECTED = 3; // Blue

    static final int defaultWidth = 200;
    static final int defaultHeight = 50;

    public StsTablePanelNew()
    {
    }

    public StsTablePanelNew( ArrayList arrayList, String[] columnNames, String[] columnTitles)
    {
        tableModel = new StsEditableTableModel( arrayList, columnNames, columnTitles);
        initializeTable();
//        initialize();
    }

    public StsTablePanelNew( ArrayList arrayList, String[] columnNames)
    {
        tableModel = new StsEditableTableModel( arrayList, columnNames, columnNames);
        initializeTable();
//        initialize();
    }

    public StsTablePanelNew( StsEditableTableModel tableModel)
    {
        this.tableModel = tableModel;
        initializeTable();
//        initialize();
    }

    public void addDeleteButtons()
    {
        tableModel.addDeleteButtons(true);
    }

    public void initializeTable()
    {
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                constructTable();
            }
        };
        StsToolkit.runWaitOnEventThread(runnable);
    }

    public void initialize()
    {
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                constructTable();
                buildPanel();
            }
        };
        StsToolkit.runWaitOnEventThread(runnable);
    }

    public void setSize(int width, int height)
    {
        this.scrollPanelWidth = width;
        this.scrollPaneHeight = height;
    }

    private void buildPanel()
    {
        selectionColor = table.getSelectionBackground();
        table.setShowGrid( true);
        table.setRequestFocusEnabled( false);
        table.getSelectionModel().addListSelectionListener( this);
        table.getModel().addTableModelListener(this);
//        table.setPreferredSize( new Dimension(scrollPanelWidth, scrollPaneHeight));
        setDeleteColumnWidth();
        JScrollPane scrollPane = new JScrollPane( table);
        scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
	    scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
 //       scrollPane.setAutoscrolls(true);
        gbc.fill = GridBagConstraints.BOTH;
        scrollPane.setPreferredSize( new Dimension(scrollPanelWidth, scrollPaneHeight));
        add(scrollPane);
//        add(table);
    }

    private void setDeleteColumnWidth()
    {
        for (int i = 0; i < tableModel.getColumnCount(); i++)
        {
            if (tableModel.getColumnName(i).equalsIgnoreCase("Delete"))
            {
                TableColumnModel tcm = table.getColumnModel();
                TableColumn col = tcm.getColumn(i);
                col.setMinWidth(50);
                col.setMaxWidth(50);
                col.setPreferredWidth(50);
                col.setResizable(false);
            }
        }
    }

    public void deleteRow()
    {
        int selectedRow = table.getSelectedRow();
        System.out.println("Deleting row " + selectedRow);
    }

    public JTable getTable()
    {
        return table;
    }

    private void constructTable()
    {
        table = new JTable(tableModel);
        tableModel.setTable(table);
    }

    public void setColumnWidth( int colIdx, int width)
    {
        TableColumn col = table.getColumnModel().getColumn(colIdx);
        col.setPreferredWidth(width);
    }

    public boolean setSelectAll( boolean selectAll)
    {
        if( !selectAll) return true;

        int nRows = tableModel.getRowCount();
        if(nRows <= 0) return false;
        int[] all = new int[nRows];
        for (int i = 0; i < all.length; i++)
           all[i] = i;
        setSelectedIndices(all);
        return true;
    }

    public void setSelectionIndex( int index)
    {
        setSelectedIndices(new int[] {index});
    }

    public void highlightRows(int[] rows, Color color)
    {
    	if(rows == null) return;
    	if(rows.length == 0) return;
    	
    	table.setSelectionBackground(color);
    	setSelectedIndices(rows);    	
    }
    
    public void resetHighlight()
    {
    	table.setSelectionBackground(selectionColor);
    	table.clearSelection();
    }

    public void setColumnSelectable( boolean val)
    {
        table.setColumnSelectionAllowed( val);
    }

    public void setRowsSelectable( boolean val)
    {
        table.setRowSelectionAllowed(val);
    }

    public void setSelectionMode( int selectionMode)
    {
        table.setSelectionMode( selectionMode);
    }

    public void setColumnSelectable()
    {
        table.setColumnSelectionAllowed(false);
    }

    public void setSelectedIndices( int[] indices)
    {
        if(indices == null) return;
        DefaultListSelectionModel model = (DefaultListSelectionModel)table.getSelectionModel();
        if(model == null) return;
        model.clearSelection();
        int nIndices = indices.length;
        model.addSelectionInterval(indices[0], indices[nIndices-1]);
    }

    public void addRow( Object row)
    {
        tableModel.addRow(row);
        repaint();
    }

    public void addRows( Object[] rowObjects)
    {
        tableModel.addRows( rowObjects);
        table.revalidate();
        repaint();
    }

    public void replaceRows( ArrayList arrayList)
    {
        if( arrayList == null || arrayList.size() == 0) return;
        tableModel.replaceRows(arrayList);
//        repaint();
    }

    public void removeRow( Object row)
    {
        int index = tableModel.removeRowObject(row);
        if(index == StsParameters.NO_MATCH) return;
//        repaint();
    }
    
    public void rowDeleteNotify(StsDeleteRowNotifyListener listener, int index)
    {
        listener.deleteRow(index);
    }

    public synchronized void addDeleteRowNotifyListener(StsDeleteRowNotifyListener listener)
    {
        if (!deleteRowNotifyListeners.contains(listener))
        {
            deleteRowNotifyListeners.add(listener);
        }
        setRowsSelectable(deleteRowNotifyListeners.size() > 0);
    }

    public synchronized void removeDeleteRowNotifyListener(StsDeleteRowNotifyListener listener)
    {
        deleteRowNotifyListeners.remove(listener);
    }
    
    public void rowSelectNotify(StsSelectRowNotifyListener listener, int[] indices)
    {
        listener.rowsSelected( indices);
    }

    public synchronized void addSelectRowNotifyListener( StsSelectRowNotifyListener listener)
    {
        if (! selectRowNotifyListeners.contains(listener))
        {
            selectRowNotifyListeners.add(listener);
        }
        setRowsSelectable( selectRowNotifyListeners.size() > 0);
    }

    public synchronized void removeSelectRowNotifyListener( StsSelectRowNotifyListener listener)
    {
        selectRowNotifyListeners.remove(listener);
        setRowsSelectable( selectRowNotifyListeners.size() > 0);
    }

    public synchronized void addListSelectionListener( ListSelectionListener l)
    {
        Vector v = listSelectionListeners == null ? new Vector(2) : (Vector)listSelectionListeners.clone();
        if(!v.contains(l))
        {
            v.addElement(l);
            listSelectionListeners = v;
        }
    }

    public synchronized void removeListSelectionListener( ListSelectionListener l)
    {
        if(listSelectionListeners != null && listSelectionListeners.contains(l))
        {
            Vector v = (Vector)listSelectionListeners.clone();
            v.removeElement(l);
            listSelectionListeners = v;
        }
    }

    public void valueChanged( ListSelectionEvent e)
    {
        table.requestFocus();
        if(!e.getValueIsAdjusting())
        {
            ListSelectionModel lsm = (ListSelectionModel)e.getSource();
            if (lsm.isSelectionEmpty())
            {
                return;
            }
            int[] indexes = table.getSelectedRows();
            Iterator itr = selectRowNotifyListeners.iterator();
            while (itr.hasNext())
            {
                StsSelectRowNotifyListener listener = (StsSelectRowNotifyListener)itr.next();
                rowSelectNotify(listener, indexes);
            }
        }
    }

    public synchronized void addTableModelListener( StsTableModelListener l)
    {
        Vector v = tableModelListeners == null ? new Vector(2) : (Vector)tableModelListeners.clone();
        if(!v.contains(l))
        {
            v.addElement(l);
            tableModelListeners = v;
        }
    }

    public synchronized void removeTableModelListener( StsTableModelListener l)
    {
        if(tableModelListeners != null && tableModelListeners.contains(l))
        {
            Vector v = (Vector)tableModelListeners.clone();
            v.removeElement(l);
            tableModelListeners = v;
        }
    }

    public void rowRemoveNotify( StsTableModelListener listener, int firstRow, int lastRow)
    {
        listener.removeRows( firstRow, lastRow);
    }

    public void tableChanged( TableModelEvent e)
    {
        int firstRow = e.getFirstRow();
        int lastRow  = e.getLastRow();

        switch (e.getType())
        {
            case TableModelEvent.DELETE:
                if( tableModelListeners != null)
                {
                    Iterator itr = tableModelListeners.iterator();
                    while (itr.hasNext())
                    {
                        StsTableModelListener listener = (StsTableModelListener)itr.next();
                        rowRemoveNotify( listener, firstRow, lastRow);
                    }
                }
                break;
        }       
    }

    public boolean isSelectedRow(int row)
    {
        for(int i = 0; i < table.getSelectedRowCount(); i++)
        {
            if(row == table.getSelectedRows()[i])
                return true;
        }
        return false;
    }

    public void setValueAt(Object value_, Object rowObject_, String columnName_)
    {
        final Object value = value_;
        final Object rowObject = rowObject_;
        final String columnName = columnName_;

        Runnable runnable = new Runnable()
        {
            public void run()
            {
                tableModel.setValueAt(value, rowObject, columnName);
                repaint();
            }
        };
        StsToolkit.runLaterOnEventThread(runnable);
    }
    
    static public void main(String[] args)
    {

        int nObjects = 3;
        ArrayList rowList = new ArrayList();
        StsSeismicBoundingBox[] rowObjects = new StsSeismicBoundingBox[nObjects];
        String[] columnNames = new String[]{"xMin", "xMax", "yMin", "yMax", "statusString"};
        String[] columnTitles = new String[]{"minX", "maxX", "minY", "maxY", "Status"};
        for(int n = 0; n < nObjects; n++)
        {
            rowObjects[n] = new StsSeismicVolume();
            rowObjects[n].initialize(n + 0.1f, n + 0.2f, n + 0.3f, n + 0.4f, n + 0.5f, n + 0.6f);
            rowObjects[n].statusString = StsSeismicBoundingBox.STATUS_OK_STR;
            rowList.add(rowObjects[n]);
        }
        StsEditableTableModel tableModel = new StsEditableTableModel(StsSeismicBoundingBox.class, columnNames, columnTitles, true);
        StsTablePanelNew tablePanel = new StsTablePanelNew(tableModel);
        tablePanel.setSize(400, 200);
        tablePanel.initialize();
        StsToolkit.createDialog(tablePanel, false);

        StsToolkit.sleep(2000);
        tablePanel.addRows( rowObjects);

        for(int n = 0; n < nObjects; n++)
        {
            StsToolkit.sleep(500);
            tablePanel.setValueAt(StsSeismicBoundingBox.STATUS_GRID_BAD_STR, rowObjects[n], "statusString");
        }
        for(int n = nObjects - 1; n >= 0; n--)
        {
            StsToolkit.sleep(500);
            tablePanel.removeRow(rowObjects[n]);
        }
        for(int n = 0; n < nObjects; n++)
        {
            StsToolkit.sleep(500);
            tablePanel.addRow(rowObjects[n]);
        }
        StsToolkit.sleep(500);
        tablePanel.addRows(rowObjects);
    }
}