package com.Sts.UI.Table;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Dec 14, 2007
 * Time: 3:39:23 PM
 * To change this template use File | Settings | File Templates.
 */

import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.table.*;
import java.lang.reflect.*;
import java.util.*;

public class StsEditableTableModel extends AbstractTableModel
{
    String[] columnTitles;
    String[] columnNames;
    ArrayList arrayList = new ArrayList();
    Field[] columnFields;
    JTable table;
    Class rowClass;
    String rowClassname;
    int nCols;
    public boolean addDeleteButtons = false;
    int nDeleteButtonColumn = -1;
    boolean[] delete;

    static final int NO_MATCH = StsParameters.NO_MATCH;

    public StsEditableTableModel()
    {
    }

    public StsEditableTableModel(ArrayList arrayList, String[] columnNames)
    {
        nCols = columnNames.length;
        this.columnNames = columnNames;
        initializeRows( arrayList);
        initializeColumnFields();
    }

    public StsEditableTableModel(ArrayList arrayList, String[] columnNames, String[] columnTitles)
    {
        nCols = columnNames.length;
        this.columnNames = columnNames;
        this.columnTitles = columnTitles;
        initializeRows(arrayList);
        initializeColumnFields();
    }

    public StsEditableTableModel(Class rowClass, String[] columnNames, String[] columnTitles)
    {
        nCols = columnNames.length;
        this.rowClass = rowClass;
        this.columnNames = columnNames;
        this.columnTitles = columnTitles;
        initializeColumnFields();
    }

    public StsEditableTableModel(Class rowClass, String[] columnNames, String[] columnTitles, boolean addDeleteButtons)
    {
        this( rowClass, columnNames, columnTitles);
        addDeleteButtons( addDeleteButtons);
    }

    public void initializeRows(ArrayList arrayList)
    {
        rowClass = arrayList.get(0).getClass();
        rowClassname = rowClass.getName();
        this.arrayList = arrayList;
    }

    public void setTable(JTable table)
    {
        this.table = table;
    }

    public boolean replaceRows( ArrayList arrayList)
    {
        if( rowClass == null || columnNames == null)
        {
            StsException.systemError( this, "replaceRows", "Can't replace rows until rowClass and columnNames initialized. ");
            return false;
        }
        if( !rowClass.isInstance( arrayList.get(0)))
        {
            StsException.systemError(this, "replaceRows", "Can't replace rows with objects of a different class. ");
            return false;
        }
        this.arrayList = arrayList;
        if( addDeleteButtons) delete = new boolean[arrayList.size()];
        return true;
    }

    public boolean addRows(Object[] rowObjects)
    {
        if(rowClass == null || columnNames == null)
        {
            StsException.systemError(this, "replaceRows", "Can't replace rows until rowClass and columnNames initialized. ");
            return false;
        }
        if( !(rowClass.isInstance( rowObjects[0])))
        {
            StsException.systemError(this, "replaceRows", "Can't replace rows with objects of a different class. ");
            return false;
        }

        for( int i = 0; i< rowObjects.length; i++)
            arrayList.add( rowObjects[i]);

        if( addDeleteButtons)
        {
            boolean[] newDelete = new boolean[rowObjects.length];
            delete = (boolean[] )StsMath.arrayAddArray(delete, newDelete);
        }
        return true;
    }

    public boolean addRow(Object rowObject)
    {
        if( !(rowClass.isInstance(rowObject)))
        {
            StsException.systemError(this, "addRow", "Cannot add instance of class " + StsToolkit.getSimpleClassname(rowObject) + " to table of " + rowClassname);
            return false;
        }
        arrayList.add( rowObject);
        if( addDeleteButtons) delete = (boolean[])StsMath.arrayAddElement(delete, false);
        return true;
    }

    public int removeRowObject(Object rowObject)
    {
        if( !(rowClass.isInstance( rowObject)))
        {
            StsException.systemError(this, "removeRowObject", "Cannot remove instance of class " + StsToolkit.getSimpleClassname(rowObject) + " to table of " + rowClassname);
            return NO_MATCH;
        }
        int index = arrayList.indexOf(rowObject);
        if(index < 0)
        {
            StsException.systemError(this, "removeRowObject", "Table doesn't contain this object");
            return NO_MATCH;
        }
        arrayList.remove(index);
        if( addDeleteButtons) delete = (boolean[] )StsMath.arrayDeleteElement(delete, index);
        return index;
    }

    private void initializeColumnFields()
    {
        columnFields = new Field[nCols];
        for(int n = 0; n < nCols; n++)
        {
            try
            {
                columnFields[n] = StsToolkit.getField(rowClass, columnNames[n]);
            }
            catch(Exception e)
            {
                StsException.systemError(this, "constructor", "Failed to find field " + columnNames[n] + " in class " + rowClass.getName());
            }
        }
    }

    public void addDeleteButtons( boolean addButton)
    {
        if( !addButton) return;
        this.addDeleteButtons = addButton;
        nDeleteButtonColumn = nCols;
        nCols++;
        columnTitles = (String[])StsMath.arrayAddElement(columnTitles, "delete");
    }

    public int getRowCount()
    {
        return arrayList.size();
    }

    public int getColumnCount()
    {
        if( columnTitles != null)
            return columnTitles.length;
        else
            return 0;
    }

    public Object getValueAt(int row, int col)
    {
        try
        {
            if(col == nDeleteButtonColumn)
                return delete[row];
            else
                return columnFields[col].get(arrayList.get(row));
        }
        catch(Exception e)
        {
            StsException.systemError(this, "getValueAt", "Failed to find value at row " + row + " col " + col);
            return null;
        }
    }

    public String getColumnName(int col)
    {
        return columnTitles[col];
    }

    public Class getColumnClass(int col)
    {
//****wrw DefaultCellRenderer gets render class here, if getValueAt is used the default render class is always a Number
//        return getValueAt(0, col).getClass();
        if(col == nDeleteButtonColumn)
            return Boolean.class;
        else
            return columnFields[col].getClass();
    }

    public boolean isCellEditable(int row, int col)
    {
        return true;
    }

    public void setValueAt(Object value, int row, int col)
    {
       try
        {
           if(col == nDeleteButtonColumn)
           {
               delete[row] = ((Boolean)value).booleanValue();
               if(delete[row]) deleteRow(row);
           }
            else
               columnFields[col].set(arrayList.get(row), value);
        }
        catch(Exception e)
        {
        }
    }

    public void setValueAt(Object aValue, Object rowObject, String columnName)
    {
       int rowIndex = -1, columnIndex = -1;
       try
       {
            rowIndex = getRowIndex(rowObject);
            if(rowIndex == -1) return;
            columnIndex = getColumnIndex(columnName);
            if(columnIndex == -1) return;
            setValueAt(aValue, rowIndex, columnIndex);
       }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "setValueAt", "Failed to set value for object " + rowClassname + "[" + rowIndex + "]." + columnName, e);
        }
    }

    private int getRowIndex(Object rowObject)
    {
        return arrayList.indexOf( rowObject);
    }

    private int getColumnIndex(String columnName)
    {
        for(int n = 0; n < nCols; n++)
            if(columnNames[n].equals(columnName)) return n;
        return -1;
    }

    private void deleteRow(int row)
    {
        fireTableRowsDeleted(row, row);
        arrayList.remove(row);
        if(addDeleteButtons) delete = (boolean[] )StsMath.arrayDeleteElement(delete, row);
        table.repaint();
    }
}
