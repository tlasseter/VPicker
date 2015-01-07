package com.Sts.UI.Table;

import com.Sts.Utilities.*;

import javax.swing.event.*;
import javax.swing.table.*;
import java.lang.reflect.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Dec 6, 2007
 * Time: 5:27:44 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsObjectTableModel implements TableModel
{
    Class rowObjectClass;
    String rowObjectClassname;
    Object[] rowObjects;
    String[] columnNames;
    String[] columnTitles = null;
    Field[] columnFields;
    int nRows;
    int nCols;
    public boolean addDeleteButtons = false;
    int nDeleteButtonColumn = -1;
    public boolean deleted = false;

    static final int NO_MATCH = StsParameters.NO_MATCH;

       /** List of listeners */
    protected EventListenerList listenerList = new EventListenerList();

    public StsObjectTableModel()
    {
    }

    public StsObjectTableModel(Object[] rowObjects, String[] columnNames)
    {
        initializeRows(rowObjects);
        initializeColumns(columnNames);
    }

    public StsObjectTableModel(Class rowObjectClass, String[] columnNames, String[] columnTitles)
    {
        this.rowObjectClass = rowObjectClass;
        this.columnTitles = columnTitles;
        initializeColumns(columnNames);
    }

    public StsObjectTableModel(Class rowObjectClass, String[] columnNames, String[] columnTitles, boolean addDeleteButtons)
    {
        this.rowObjectClass = rowObjectClass;
        this.columnTitles = columnTitles;
        this.addDeleteButtons = addDeleteButtons;
        initializeColumns(columnNames);
    }

    public void initializeRows(Object[] rowObjects)
    {
        rowObjectClass = rowObjects[0].getClass();
        rowObjectClassname = rowObjectClass.getName();
        setRowObjects(rowObjects);
    }

    private void setRowObjects(Object[] rowObjects)
    {
        if(rowObjects == null)
        {
            this.rowObjects = rowObjects;
            nRows = 0;
        }
        else
        {
            this.rowObjects = rowObjects;
            nRows = rowObjects.length;
        }
    }

    public boolean replaceRows(Object[] rowObjects)
    {
        if(rowObjectClass == null || columnNames == null)
        {
            StsException.systemError(this, "replaceRows", "Can't replace rows until rowObjectClass and columnNames initialized. ");
            return false;
        }
        if(rowObjects[0].getClass() != rowObjectClass)
        {
            StsException.systemError(this, "replaceRows", "Can't replace rows with objects of a different class. ");
            return false;
        }
        setRowObjects(rowObjects);
        return true;
    }

    public boolean addRows(Object[] rowObjects)
    {
        if(rowObjectClass == null || columnNames == null)
        {
            StsException.systemError(this, "replaceRows", "Can't replace rows until rowObjectClass and columnNames initialized. ");
            return false;
        }
        if(rowObjects[0].getClass() != rowObjectClass)
        {
            StsException.systemError(this, "replaceRows", "Can't replace rows with objects of a different class. ");
            return false;
        }
        rowObjects = (Object[]) StsMath.arrayAddArray(this.rowObjects, rowObjects);
        setRowObjects(rowObjects);
        return true;
    }

    public boolean addRow(Object rowObject)
    {
        if(rowObject.getClass() != rowObjectClass)
        {
            StsException.systemError(this, "addRow", "Cannot add instance of class " + StsToolkit.getSimpleClassname(rowObject) + " to table of " + rowObjectClassname);
            return false;
        }
        setRowObjects((Object[])StsMath.arrayAddElement(rowObjects, rowObject));
        return true;
    }

    public void removeAllRows()
    {
        setRowObjects(null);
    }

    public int removeRowObject(Object rowObject)
    {
        if(rowObject.getClass() != rowObjectClass)
        {
            StsException.systemError(this, "addRow", "Cannot add instance of class " + StsToolkit.getSimpleClassname(rowObject) + " to table of " + rowObjectClassname);
            return NO_MATCH;
        }
        int index = StsMath.index(rowObjects, rowObject);
        if(index == NO_MATCH)
        {
            StsException.systemError(this, "removeRowObject", "Table doesn't contain this object");
            return NO_MATCH;
        }
        setRowObjects((Object[])StsMath.arrayDeleteElement(rowObjects, index));
        return index;
    }

    public void initializeColumns(String[] columnNames)
    {
        this.columnNames = columnNames;
        nCols = columnNames.length;
        columnFields = new Field[nCols];
        for(int n = 0; n < nCols; n++)
        {
            try
            {
                columnFields[n] = rowObjectClass.getDeclaredField(columnNames[n]);
            }
            catch(Exception e)
            {
                StsException.systemError(this, "constructor", "Failed to find field " + columnNames[n] + " in class " + rowObjectClass.getName());
            }
        }
        if(addDeleteButtons) addDeleteButton();
    }

    public void addDeleteButton()
    {
        nDeleteButtonColumn = nCols;
        nCols++;
        columnNames = (String[])StsMath.arrayAddElement(columnNames, "Delete");
        try
        {
            Field deleteField = StsObjectTableModel.class.getDeclaredField("deleted");
            columnFields = (Field[])StsMath.arrayAddElement(columnFields, deleteField);
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "addDeleteButtonColumn", e);
        }
    }

    /**
     * Returns the number of rows in the model. A
     * <code>JTable</code> uses this method to determine how many rows it
     * should display.  This method should be quick, as it
     * is called frequently during rendering.
     *
     * @return the number of rows in the model
     * @see #getColumnCount
     */
    public int getRowCount()
    {
        return nRows;
    }

    /**
     * Returns the number of columns in the model. A
     * <code>JTable</code> uses this method to determine how many columns it
     * should create and display by default.
     *
     * @return the number of columns in the model
     * @see #getRowCount
     */
    public int getColumnCount()
    {
        return nCols;
    }

    /**
     * Returns the name of the column at <code>columnIndex</code>.  This is used
     * to initialize the table's column header name.  Note: this name does
     * not need to be unique; two columns in a table can have the same name.
     *
     * @param	columnIndex	the index of the column
     * @return  the name of the column
     */
    public String getColumnName(int columnIndex)
    {
        if( columnTitles != null && columnIndex < columnTitles.length)
            return columnTitles[columnIndex];
        else
            return columnNames[columnIndex];
    }

    /**
     * Returns the most specific superclass for all the cell values
     * in the column.  This is used by the <code>JTable</code> to set up a
     * default renderer and editor for the column.
     *
     * @param columnIndex  the index of the column
     * @return the common ancestor class of the object values in the model.
     */
    public Class<?> getColumnClass(int columnIndex)
    {
        return columnFields[columnIndex].getClass();
    }

    /**
     * Returns true if the cell at <code>rowIndex</code> and
     * <code>columnIndex</code>
     * is editable.  Otherwise, <code>setValueAt</code> on the cell will not
     * change the value of that cell.
     *
     * @param	rowIndex	the row whose value to be queried
     * @param	columnIndex	the column whose value to be queried
     * @return	true if the cell is editable
     * @see #setValueAt
     */
    public boolean isCellEditable(int rowIndex, int columnIndex)
    {
        return columnIndex == nDeleteButtonColumn;
    }

    /**
     * Returns the value for the cell at <code>columnIndex</code> and
     * <code>rowIndex</code>.
     *
     * @param	rowIndex	the row whose value is to be queried
     * @param	columnIndex 	the column whose value is to be queried
     * @return	the value Object at the specified cell
     */
    public Object getValueAt(int rowIndex, int columnIndex)
    {
        try
        {
            if(columnIndex == nDeleteButtonColumn) return null;
            return columnFields[columnIndex].get(rowObjects[rowIndex]);
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "getValueAt", "Failed to get value for object " + rowObjectClassname + "[" + rowIndex + "]." + columnNames[columnIndex], e);
            return null;
        }
    }

    /**
     * Sets the value in the cell at <code>columnIndex</code> and
     * <code>rowIndex</code> to <code>aValue</code>.
     *
     * @param	aValue		 the new value
     * @param	rowIndex	 the row whose value is to be changed
     * @param	columnIndex 	 the column whose value is to be changed
     * @see #getValueAt
     * @see #isCellEditable
     */
    public void setValueAt(Object aValue, int rowIndex, int columnIndex)
    {
       try
       {
           if(columnIndex == nDeleteButtonColumn) return;
           columnFields[columnIndex].set(rowObjects[rowIndex], aValue);
       }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "setValueAt", "Failed to set value for object " + rowObjectClassname + "[" + rowIndex + "]." + columnNames[columnIndex], e);
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
            StsException.outputWarningException(this, "setValueAt", "Failed to set value for object " + rowObjectClassname + "[" + rowIndex + "]." + columnName, e);
        }
    }

    private int getRowIndex(Object rowObject)
    {
        for(int n = 0; n < nRows; n++)
            if(rowObjects[n] == rowObject) return n;
        return -1;
    }


    private int getColumnIndex(String columnName)
    {
        for(int n = 0; n < nCols; n++)
            if(columnNames[n].equals(columnName)) return n;
        return -1;
    }

    /**
     * Adds a listener to the list that is notified each time a change
     * to the data model occurs.
     *
     * @param	l		the TableModelListener
     */
    public void addTableModelListener(TableModelListener l)
    {
	    listenerList.add(TableModelListener.class, l);
    }

    /**
     * Removes a listener from the list that's notified each time a
     * change to the data model occurs.
     *
     * @param	l		the TableModelListener
     */
    public void removeTableModelListener(TableModelListener l)
    {
	    listenerList.remove(TableModelListener.class, l);
    }
}
