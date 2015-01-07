package com.Sts.Actions.Wizards.TraceDefinition;

import com.Sts.Types.*;
import com.Sts.Utilities.*;

import javax.swing.table.*;
import java.text.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2006</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */
public class StsSegyTraceTableModel extends AbstractTableModel
{
	private Object[] data = null;
	private int startTrace;
	private int traceColumnCount;
	private int xyScalePos = 0;
    private int xyScaleFmt = 0;
    private int edScalePos = 0;
    private int edScaleFmt = 0;
	private byte[][] traceBinaryHeader = null;
	private boolean littleEndian;

	private DecimalFormat labelFormat = new DecimalFormat("####.#");
	private DecimalFormat intFormat = new DecimalFormat("#");
	private float scalar = 1.0f;
	private int fieldPos = 0;
	private int format;

	public StsSegyTraceTableModel(Object[] data, int startTrace, int traceColumnCount, int xyScalePos, int xyScaleFmt,
                                     int edScalePos, int edScaleFmt, byte[][] traceBinaryHeader, boolean littleEndian)
	{
		super();
		this.data = data;
		this.startTrace = startTrace;
		this.traceColumnCount = traceColumnCount;
		this.xyScalePos = xyScalePos;
		this.xyScaleFmt = xyScaleFmt;
        this.edScalePos = edScalePos;
		this.edScaleFmt = edScaleFmt;
		this.traceBinaryHeader = traceBinaryHeader;
		this.littleEndian = littleEndian;
	}

	/**
	 * Returns the number of columns in the model.
	 *
	 * @return the number of columns in the model
	 */
	public int getColumnCount()
	{
		return traceColumnCount + 1;
	}

	public String getColumnName(int columnIndex)
	{
		if (columnIndex == 0)
		{
			return "Attribute";
		}
		return "T" + (startTrace + columnIndex - 1);
	}

	public boolean isCellEditable(int rowIndex, int columnIndex)
	{
		return false;
	}

	public Class getColumnClass(int columnIndex)
	{
    	return String.class;
	}

	/**
	 * Returns the number of rows in the model.
	 *
	 * @return the number of rows in the model
	 */
	public int getRowCount()
	{
		return data.length;
	}

	/**
	 * Returns the value for the cell at <code>columnIndex</code> and <code>rowIndex</code>.
	 *
	 * @param rowIndex the row whose value is to be queried
	 * @param columnIndex the column whose value is to be queried
	 * @return the value Object at the specified cell
	 */
	public Object getValueAt(int rowIndex, int columnIndex)
	{
		StsSEGYFormatRec record = (StsSEGYFormatRec)data[rowIndex];
        int scalePos = -1;
        int scaleFmt = 0;
	    if (columnIndex == 0)
		{
			return record.getUserName();
		}

		scalar = 1.0f;
		fieldPos = record.getLoc();
		format = record.getFormat();
        if(record.getApplyScalar().equals("CO-SCAL"))
        {
            scalePos = xyScalePos;
            scaleFmt = xyScaleFmt;
        }
        else if(record.getApplyScalar().equals("ED-SCAL"))
        {
            scalePos = edScalePos;
            scaleFmt = edScaleFmt;
        }
        try
        {
            if(scalePos != -1)
            {
                switch(scaleFmt)
                {
                    case StsSEGYFormat.IBMFLT: // IBM Float
                        scalar = StsMath.convertIBMFloatBytes(traceBinaryHeader[columnIndex - 1], scalePos, littleEndian);
                        break;
                    case StsSEGYFormat.IEEEFLT: // IEEE Float
                        scalar = Float.intBitsToFloat(StsMath.convertIntBytes(traceBinaryHeader[columnIndex - 1], scalePos, littleEndian));
                        break;
                    case StsSEGYFormat.INT4: // Integer 4
                        scalar = (float)StsMath.convertIntBytes(traceBinaryHeader[columnIndex - 1], scalePos, littleEndian);
                        break;
                    case StsSEGYFormat.INT2: // Integer 2
                        scalar = (float)StsMath.convertBytesToShort(traceBinaryHeader[columnIndex - 1], scalePos, littleEndian);
                        break;
                    default:
                        scalar = 1.0f;
                }
                if (scalar < 0.0f) scalar = -1/scalar;
                else if(scalar == 0.0f) scalar = 1.0f;
            }
            switch (format)
            {
                case StsSEGYFormat.IBMFLT: // IBM Float
                    return labelFormat.format(StsMath.convertIBMFloatBytes(traceBinaryHeader[columnIndex - 1], fieldPos, littleEndian) * scalar);
                case StsSEGYFormat.IEEEFLT: // IEEE Float
                    return labelFormat.format(Float.intBitsToFloat(StsMath.convertIntBytes(traceBinaryHeader[columnIndex - 1], fieldPos, littleEndian)) * scalar);
                case StsSEGYFormat.INT4: // Integer 4
                    return intFormat.format(StsMath.convertIntBytes(traceBinaryHeader[columnIndex - 1], fieldPos, littleEndian) * scalar);
                case StsSEGYFormat.INT2: // Integer 2
                    return intFormat.format(StsMath.convertBytesToShort(traceBinaryHeader[columnIndex - 1], fieldPos, littleEndian) * scalar);
                default:
                    return "----";
            }
        }
        catch(Exception e)
        {
            StsException.systemError(this, "getValueAt", "columnIndex: " + columnIndex + " scalePos: " + scalePos + " scalar: " + scalar + " fieldPos" + fieldPos) ;
            return "----";    
        }

	}
}
