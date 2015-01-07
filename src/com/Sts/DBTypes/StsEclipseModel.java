package com.Sts.DBTypes;

import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.StsCursor;
import com.Sts.Utilities.*;
import com.Sts.Interfaces.StsTreeObjectI;

import javax.swing.*;
import java.text.*;
import java.util.*;
import java.io.*;
import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Feb 13, 2009
 * Time: 12:00:21 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsEclipseModel extends StsMainTimeObject implements StsTreeObjectI
{
    public StsBuiltModel builtModel;
    public String pathname;
    public StsPropertyType porosity;
    public StsPropertyType permX;
    public StsPropertyType permY;
    public StsPropertyType permZ;
    public StsPropertyType actnum ;
    public StsPropertyType depth;
    public StsPropertyType poreVolume;
    public StsPropertyType tranX;
    public StsPropertyType tranY;
    public StsPropertyType tranZ;

    public StsDynamicPropertyType pressure;
    public StsDynamicPropertyType wsat;
        
    public String[] dates;

    transient public String inputPathname;
    transient public String outputPathname;

    transient public StsPropertyType none = StsPropertyType.constructor(StsPropertyType.none, false);
    transient public StsPropertyType layerColorType = StsPropertyType.constructor(StsPropertyType.layerColor, false);
    transient public StsPropertyType indexMap = StsPropertyType.constructor(StsPropertyType.indexMap, false);

    transient StsPropertyType[] staticPropertyTypes;
    transient StsDynamicPropertyType[] dynamicPropertyTypes;
    transient StsPropertyType[] displayPropertyTypes;

    transient private ArrayList<StsPropertyVolume> blockPropertyVolumes = new ArrayList<StsPropertyVolume>();
    transient private ArrayList<StsDynamicPropertyVolume> dynamicBlockPropertyVolumes = new ArrayList<StsDynamicPropertyVolume>();

    transient private StsPropertyType displayPropertyType = StsPropertyType.propertyTypeNone;
    transient public int timeSteps = 0;
    transient boolean attributeChanged = true;

    transient DefaultBoundedRangeModel rangeModel = new DefaultBoundedRangeModel(0, 0, 0, 83);
    transient int currentDateIndex = -1;
    transient String currentDate = null;

    static StsBooleanFieldBean visibleBean = new StsBooleanFieldBean(StsEclipseModel.class, "isVisible", "Enable");
    static StsComboBoxFieldBean propertyTypeBean = new StsComboBoxFieldBean(StsEclipseModel.class, "displayPropertyType", "Display", "displayPropertyTypes");
    static StsButtonFieldBean movieBean = new StsButtonFieldBean("Movie", "Play dynamic data selection as a movie.", StsEclipseModel.class, "playMovie");
    static StsEditableColorscaleFieldBean colorscaleBean = new StsEditableColorscaleFieldBean(StsEclipseModel.class, "colorscale");
    static StsComboBoxFieldBean currentDateBean = new StsComboBoxFieldBean(StsEclipseModel.class, "currentDate", "Display date", "dates");
    //static StsSliderFieldBean dateSliderBean = new StsSliderFieldBean(StsEclipseModel.class, "rangeModel", "rangeChanged");

    static public SimpleDateFormat dateTimeFormat;
    static public StsFieldBean[] displayFields = null;
    static protected StsObjectPanel objectPanel = null;

    static Thread movieThread = null;

    public StsEclipseModel()
    {
    }

    public StsEclipseModel(StsBuiltModel builtModel)
    {
        super(false);
        this.builtModel = builtModel;
        dateTimeFormat = new SimpleDateFormat(currentModel.getProject().getTimeDateFormatString());
        addToModel();
        setName("EclipseModel-" + getIndex());
        createFileDirectories(builtModel);
        constructPropertyTypes();
        initializePropertyTypes();
    }

    private void createFileDirectories(StsBuiltModel builtModel)
    {
        pathname = builtModel.getFilesPathname() + name + File.separator;
        dbFieldChanged("pathname", pathname);
        initializeFileDirectories();
        checkClearDirectories();
    }

    private void initializeFileDirectories()
    {
        inputPathname = pathname + "Input" + File.separator ;
        outputPathname = pathname + "Output" + File.separator ;
    }

    private void checkClearDirectories()
    {
        StsFile.clearDirectoryAndFiles(pathname);
        StsFile.checkCreateDirectory(inputPathname);
        StsFile.checkCreateDirectory(outputPathname);
    }

    public boolean deleteModel(StsModel model)
    {
		StsActionManager actionManager = null;

        StsCursor cursor = new StsCursor(model.win3d, Cursor.WAIT_CURSOR);
        try
        {
			actionManager = model.mainWindowActionManager;

            StsBlock[] blocks = builtModel.getBlocks();
            for(StsBlock block : blocks)
                block.clearTransientArrays();

            StsZone[] zones = (StsZone[])model.getCastObjectList(StsZone.class);
            for(StsZone zone : zones)
                zone.clearTransientArrays();

            delete();
            model.refreshObjectPanel();
            return true;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "deleteModel", e);
            return false;
        }
        finally
        {
             cursor.restoreCursor();
			 if(actionManager != null) model.getDatabase().saveTransactions();
             model.enableDisplay();
        }
    }

    public boolean initialize(StsModel model)
    {
        dateTimeFormat = new SimpleDateFormat(model.getProject().getTimeDateFormatString());
        initializeFileDirectories();
        initializePropertyTypes();
        if(dates != null)
        {
            currentDateIndex = 0;
            currentDate = dates[0];
            currentDateBean.setValueObject(currentDate);
        }
        return true;
    }

    private void constructPropertyTypes()
    {
        porosity = StsPropertyType.constructor(StsPropertyType.porosity);
        permX = StsPropertyType.constructor(StsPropertyType.permX);
        permY = StsPropertyType.constructor(StsPropertyType.permY);
        permZ = StsPropertyType.constructor(StsPropertyType.permZ);
        actnum = StsPropertyType.constructor(StsPropertyType.actnum);
        poreVolume = StsPropertyType.constructor(StsPropertyType.poreVolume);
        depth = StsPropertyType.constructor(StsPropertyType.depth);
        tranX = StsPropertyType.constructor(StsPropertyType.tranX);
        tranY = StsPropertyType.constructor(StsPropertyType.tranY);
        tranZ = StsPropertyType.constructor(StsPropertyType.tranZ);
        pressure = StsDynamicPropertyType.constructor(StsPropertyType.pressure);
        wsat = StsDynamicPropertyType.constructor(StsPropertyType.wsat);
    }

    public void initializePropertyTypes()
    {
        displayPropertyType = none;
        staticPropertyTypes = new StsPropertyType[] { porosity, permX, permY, permZ, poreVolume, actnum, depth, tranX, tranY, tranZ, indexMap };
        displayPropertyTypes = new StsPropertyType[] { none, porosity, permX, permY, permZ, poreVolume, actnum, depth, tranX, tranY, tranZ, pressure, wsat };
        dynamicPropertyTypes = new StsDynamicPropertyType[] { pressure, wsat };
    }

    public StsPropertyVolume createBlockPropertyVolume(StsBlock block, StsPropertyType propertyType, byte distributionType)
    {
        StsPropertyVolume blockPropertyVolume = new StsPropertyVolume(propertyType, block, distributionType, outputPathname);
        addBlockPropertyVolume(blockPropertyVolume);
        return blockPropertyVolume;
    }

    public void addBlockPropertyVolume(StsPropertyVolume blockPropertyVolume)
    {
        blockPropertyVolumes.add(blockPropertyVolume);
    }

    public StsPropertyVolume getBlockPropertyVolume(StsBlock block, StsPropertyType propertyType)
    {
        for(StsPropertyVolume volume : blockPropertyVolumes)
            if(volume.matches(block, propertyType)) return volume;
        return null;
    }

    public StsPropertyVolume getLoadBlockPropertyVolume(StsBlock block, StsPropertyType propertyType)
    {
        if(propertyType instanceof StsDynamicPropertyType)
        {
            if(currentDate == null) return null;
            // if(currentDateIndex == -1) return null;
            // String currentDate = dates[currentDateIndex];
            for(StsDynamicPropertyVolume volume : dynamicBlockPropertyVolumes)
                if(volume.matches(block, propertyType, currentDate)) return volume;
            StsDynamicPropertyVolume propertyVolume = new StsDynamicPropertyVolume((StsDynamicPropertyType)propertyType, outputPathname, block, currentDate);
            if(!propertyVolume.loadDataFile())
            {
                StsException.systemError(this, "getLoadBlockPropertyVolume", "Failed to find dynamic file: " + propertyVolume.filename);
                return null;
            }
            dynamicBlockPropertyVolumes.add(propertyVolume);
            return propertyVolume;
        }
        else if(propertyType instanceof StsPropertyType)
        {
            for(StsPropertyVolume volume : blockPropertyVolumes)
                if(volume.matches(block, propertyType)) return volume;
            StsPropertyVolume propertyVolume = new StsPropertyVolume(propertyType, block, outputPathname);
            if(!propertyVolume.loadDataFile())
            {
                StsException.systemError(this, "getLoadBlockPropertyVolume", "Failed to find static file: " + propertyVolume.filename);
                return null;
            }
            blockPropertyVolumes.add(propertyVolume);
            propertyVolume.loadDataFile();
            return propertyVolume;
        }
        else
            return null;
    }

    public void adjustRange(StsPropertyType propertyType, StsBlock block)
    {
        StsPropertyVolume propertyVolume = getBlockPropertyVolume(block, propertyType);
        propertyType.adjustRange(propertyVolume.valueMin,  propertyVolume.valueMax);
    }

    public void outputDataFile(StsPropertyType propertyType, StsBlock block)
    {
        StsPropertyVolume propertyVolume = getBlockPropertyVolume(block, propertyType);
        StsFloatVector floatVector = propertyVolume.computeFloatVector();
        floatVector.setMinValue(propertyType.getValueMin());
        floatVector.setMaxValue(propertyType.getValueMax());
        propertyVolume.outputDataFile(floatVector);
    }

    public String getOutputPathname()
    {
        return builtModel.getFilesPathname() + name + File.separator + "Output";
    }

    public String getInputPathname()
    {
        return builtModel.getFilesPathname() + name + File.separator + "Input";
    }

    public StsPropertyType getDisplayPropertyType()
    {
        return displayPropertyType;
    }

    public void setDisplayPropertyType(StsPropertyType displayPropertyType)
    {
        if(this.displayPropertyType == displayPropertyType) return;
        this.displayPropertyType = displayPropertyType;
        loadPropertyVolumes();
        colorscaleBean.setValueObject(displayPropertyType.colorscale);
    }

    public StsColorscale getColorscale()
    {
        if(displayPropertyType == null) return null;
        return displayPropertyType.colorscale;
    }

    public void setColorscale(StsColorscale colorscale)
    {
        if(displayPropertyType == null) return;
        displayPropertyType.colorscale = colorscale;    
    }

    private void loadPropertyVolumes()
    {
        StsBlock blocks[] = builtModel.getBlocks();
        for(StsBlock block : blocks)
        {
            StsPropertyVolume blockPropertyVolume = getLoadBlockPropertyVolume(block, displayPropertyType);
            block.setCurrentPropertyVolume(blockPropertyVolume);
        }
        currentModel.viewObjectRepaint(this, displayPropertyType);
    }

    public ArrayList<StsPropertyVolume> getBlockPropertyVolumes()
    {
        return blockPropertyVolumes;
    }

    static public void playMovie()
    {
        if(movieThread != null) movieThread.interrupt();
        movieThread = StsToolkit.checkRunRunnable
        (
            new Runnable()
            {
                public void run()
                {
                    StsEclipseModel eclipseModel = (StsEclipseModel)getCurrentModel().getCurrentObject(StsEclipseModel.class);
                    if(eclipseModel == null) return;
                    eclipseModel.doPlayMovie();
                }
            }
        );
    }

    public void doPlayMovie()
    {
        if(!(displayPropertyType instanceof StsDynamicPropertyType)) return;
        if(dates == null) return;
        int nDates = dates.length;
        for(int n = 0; n < nDates; n++)
        {
            final String nextDate = dates[n];
            StsToolkit.runWaitOnEventThread
            (
                new Runnable()
                {
                    public void run()
                    {
                        currentDateBean.setValueObject(nextDate);
                    }
                }
            );
            StsToolkit.sleep(2000);
        }
    }

    public boolean delete()
    {
        super.delete();
        if(pathname != null) StsFile.clearDirectoryAndFiles(pathname);
        return true;
    }

    public boolean readFile(StsAbstractFile restartFile)
    {
        File outputDirectory = new File(outputPathname);

        try
        {
            if(!outputDirectory.exists())
                outputDirectory.mkdirs();
        }
        catch(Exception e)
        {
            StsException.systemError(this, "outputDataFiles", "Failed to create directory " + outputPathname);
        }

        for(StsDynamicPropertyType propertyType : dynamicPropertyTypes)
        {
            String eclipseName = propertyType.eclipseName;
            File attributeDirectory = new File(outputPathname + eclipseName);
            attributeDirectory.mkdirs();
        }
        long offset = 0;
        byte[] temp = new byte[4];
        EclipseRecord record;
        timeSteps = 0;
        String date = "0_0_0";
        InputStream is = null;
        try
        {
            is = restartFile.getInputStream();
            BufferedInputStream bis = new BufferedInputStream(is);
            DataInputStream dis = new DataInputStream(bis);
            ArrayList<String> dateVector = new ArrayList<String>();
            // Process the file.
            int nRecord = 1;
            while(true)
            {
                record = eclRead(dis, offset, nRecord);
                if(record == null) break;
                nRecord++;
                if(record.header.contains("INTEHEAD"))
                {
                    int day = record.getIntegerAt(64);
                    int month = record.getIntegerAt(65);
                    int year = record.getIntegerAt(66);
                    try
                    {
                        Calendar calendar = Calendar.getInstance();
                        calendar.set(year, month, day, 0, 0, 0);
                        date = year + "_" + month + "_" + day;
                        dateVector.add(date);
                    }
                    catch(Exception e)
                    {
                        StsException.systemError(this, "readFile", "Failed to parse date: " + date);
                    }
                }
                else
                {
                    for(StsDynamicPropertyType propertyType : dynamicPropertyTypes)
                    {
                        if(record.header.contains(propertyType.eclipseName))
                            outputRecord(record, propertyType, date);
                    }
                }
            }
            int nTimes = dateVector.size();
            dates = new String[nTimes];
            dateVector.toArray(dates);
            dbFieldChanged("dates", dates);
            date = dates[0];
            
            for(StsDynamicPropertyType propertyType : dynamicPropertyTypes)
                propertyType.rangeChanged();
        }
        catch (IOException ex)
        {
            System.out.println("Exception processing Eclipse file: " + ex);
            return false;
        }
        finally
        {
            try
            {
                if(is != null) is.close();
            }
            catch(Exception e)
            {
            }
        }
        return true;
    }

    private boolean outputRecord(EclipseRecord record, StsDynamicPropertyType propertyType, String date)
    {
        byte[] bytes = record.getData();
        int nFloats = bytes.length/4;
        float[] floats = new float[nFloats];
        StsMath.convertIEEEBytesToFloats(bytes, floats, nFloats, false);
        outputDynamicDataFile(propertyType, date, floats);
        return true;
    }

    private void outputDynamicDataFile(StsDynamicPropertyType propertyType, String date, float[] values)
    {
        StsBlock blocks[] = builtModel.getBlocks();
        for(StsBlock block : blocks)
        {
            StsDynamicPropertyVolume propertyVolume = new StsDynamicPropertyVolume(propertyType, outputPathname, block, StsPropertyVolume.PROP_VARIABLE_XYZ, date);
            float[] indexMapValues = this.getLoadBlockPropertyVolume(block, indexMap).getValues();
            StsFloatVector floatVector = propertyVolume.computeDynamicFloatVector(values, indexMapValues);
            propertyType.adjustRange(floatVector.getMinValue(), floatVector.getMaxValue());
            propertyVolume.outputDataFile(floatVector);
        }
    }

    private EclipseRecord eclRead(DataInputStream dis, long offset, int nRecord)
    {
        try
        {
            // Read Header Block
            EclipseRecord record = new EclipseRecord();
            record.setOffset(offset);
            // UnFormatted file
            byte[] keyArray = new byte[ECL_STRING_LENGTH];
            byte[] fmtArray = new byte[ECL_TYPE_LENGTH];

            String hdr = null;
            String fmt;
            int recSize = dis.readInt();    // Start Header Block
            int read = 0;
            while(!isValidKey(hdr))
            {
                read += dis.read(keyArray, 0, ECL_STRING_LENGTH);
                hdr = new String(keyArray);
                if(read == recSize)
                {
                    dis.readInt();
                    System.out.println("Encountered 8 byte header record (" + hdr + ")");
                    return record;
                }
            }
            StsException.systemDebug(this, "eclRead", "Record  " + nRecord + " is " + hdr);
            record.setHeader(hdr);
            offset += ECL_STRING_LENGTH;

            int size = dis.readInt();
            record.setSize(size);
            offset += 4;
            read += 4;
            if(read == recSize)
            {
                dis.readInt();
                return record;
            }

            read += dis.read(fmtArray, 0, ECL_TYPE_LENGTH);
            fmt = new String(fmtArray);
            record.setEclFmt(record.getEclFmtFromString(fmt));
            offset += ECL_TYPE_LENGTH;

            recSize = dis.readInt(); // End Header Block
            if(size == 0)
            {
                System.out.println("Encountered 0 byte data block for (" + hdr + ")");
                return record;
            }
            // Read Data Block
            int blockSize = record.getBlockSize();
            read = 0;
            if((record.getEclFmt() == ECL_CHAR_TYPE) || (record.getEclFmt() == ECL_MESS_TYPE))
            {
                recSize = dis.readInt(); // Start Data Block
                int blocks = size /  blockSize + (size % blockSize == 0 ? 0 : 1);
                for (int ib = 0; ib < blocks; ib++)
                {
                    int read_elm = StsMath.min2((ib + 1) * blockSize , record.getSize() - ib * blockSize);
                    int ir;
                    for (ir = 0; ir < read_elm; ir++)
                    {
                        read += dis.read(keyArray, 0, ECL_STRING_LENGTH);
                        fmt = new String(keyArray);
                        // Accumulate the blocks in a single string
                        //fread(&ecl_kw->data[(ib * blocksize + ir) * ecl_kw->sizeof_ctype] , 1 , ECL_STRING_LENGTH , stream);
                    }
                }
                recSize = dis.readInt(); // End Data Block
            }
            else
            {
                int bufSize = record.getSize() * record.getTypeSize();
                int maxBlockSize = BLOCKSIZE_NUMERIC * record.getTypeSize();
                int blocks = 1;
                int remain = 0;
                if(record.getSize() > BLOCKSIZE_NUMERIC)
                {
                    remain = bufSize%maxBlockSize;
                    blocks = bufSize/maxBlockSize + 1;
                }
                byte[] buffer = new byte[bufSize];
                for(int ib=0; ib<blocks; ib++)
                {
                    bufSize = dis.readInt();
                    read += dis.read(buffer, read, bufSize);
                    bufSize = dis.readInt();
                    //if(read == (record.getSize() * record.getTypeSize()))
                    //    break;
                }
                record.setData(buffer);
            }
            return record;
        }
        catch(IOException ex)
        {
            if(ex instanceof EOFException)
            {
                System.out.println("Finished processing file.....");
                return null;
            }
            System.out.println("Exception processing Eclipse file header: " + ex);
            return null;
        }
    }

    public String getBornTimeString() { return dateTimeFormat.format(bornDate); }
    public String getDeathTimeString() { return dateTimeFormat.format(deathDate); }

    public String getDateTimeString(long dateTime)
    {
        return dateTimeFormat.format(dateTime);
    }
    public boolean anyDependencies() { return false; }

    public void treeObjectSelected()
    {
        currentModel.getCreateStsClass("com.Sts.DBTypes.StsEclipseModel").selected(this);
    }

     public StsFieldBean[] getDisplayFields()
     {
         if(displayFields == null)
         {
            displayFields = new StsFieldBean[]
            {
                 visibleBean,
                 propertyTypeBean,
                 currentDateBean,
                 movieBean,
                 colorscaleBean
             };
         }
         // dateSliderBean.setRangeModel(rangeModel);
         // dateSliderBean.setEditable(true);
         return displayFields;
     }
/*
     public StsFieldBean[] getDisplayFields()
     {
         if(dates != null)
         {
            int nTimes = dates.length;
            rangeModel.setRangeProperties(0, 0, 0, nTimes-1, false);
            dateSliderBean.setRangeModel(rangeModel);
            if(currentDateIndex == -1)
                currentDateIndex = 0;
            displayFields = new StsFieldBean[]
            {
                 visibleBean,
                 propertyTypeBean,
                 currentDateBean,
                 dateSliderBean,
                 colorscaleBean
             };
         }

         else if(displayFields == null)
         {
             colorscaleBean = new StsEditableColorscaleFieldBean(StsEclipseModel.class, "colorscale");
             displayFields = new StsFieldBean[]
             {
                 visibleBean,
                 propertyTypeBean,
                 colorscaleBean
             };
         }
         return displayFields;
     }
*/
    /*
    public void rangeChanged()
    {
        int currentDateIndex = rangeModel.getValue();
        currentDateBean.setValue(dates[currentDateIndex]);
        if(rangeModel.getValueIsAdjusting()) return;
        if(displayPropertyType instanceof StsDynamicPropertyType)
            loadPropertyVolumes();
    }
    */
    public String getCurrentDate()
    {
        if(dates == null) return "none";
        return currentDate;
        // return dates[currentDateIndex];
    }

    public void setCurrentDate(String date)
    {
        if(dates == null) return;
        if(currentDate == date) return;
        currentDate = date;
        if(displayPropertyType instanceof StsDynamicPropertyType)
            loadPropertyVolumes();
//        currentDateIndex = 0;
    }

    public String[] getDates() { return dates; }

     private Object[] getDisplayPropertyTypes()
     {
        return displayPropertyTypes;
     }

    public StsFieldBean[] getPropertyFields()
    {
         return null;
    }

    public Object[] getChildren() { return new Object[0]; }
    public StsObjectPanel getObjectPanel()
    {
        if(objectPanel == null) objectPanel = StsObjectPanel.constructor(this, true);
        return objectPanel;
    }
/*
    public long[] getTimes()
    {
        return times;
    }
*/
    public StsFloatVector getVectorNamed(String attName, EclipseRecord record)
    {
        if(record.getHeader().equalsIgnoreCase(attName))
        {
            StsFloatVector vec = record.getFloatVector();
            vec.setMinMax();
            return vec;
        }
        else
        {
            System.out.println("Not a " + attName + " record....");
            return null;
        }
    }

    static final String ECL_TYPE_NAME_CHAR =   "CHAR";
    static final String ECL_TYPE_NAME_FLOAT = "REAL";
    static final String  ECL_TYPE_NAME_INT = "INTE";
    static final String  ECL_TYPE_NAME_DOUBLE = "DOUB";
    static final String  ECL_TYPE_NAME_BOOL = "LOGI";
    static final String  ECL_TYPE_NAME_MESSAGE = "MESS";
    static final String[] ECL_TYPE_NAMES = {ECL_TYPE_NAME_CHAR, ECL_TYPE_NAME_FLOAT, ECL_TYPE_NAME_INT, ECL_TYPE_NAME_DOUBLE,
                                            ECL_TYPE_NAME_BOOL, ECL_TYPE_NAME_MESSAGE};

    static public final byte ECL_CHAR_TYPE = 1;
    static public final byte ECL_INT_TYPE = 2;
    static public final byte ECL_FLOAT_TYPE = 3;
    static public final byte ECL_DOUBLE_TYPE = 4;
    static public final byte ECL_BOOL_TYPE = 5;
    static public final byte ECL_MESS_TYPE = 6;
    static public final byte[] ECL_TYPES = {ECL_CHAR_TYPE, ECL_FLOAT_TYPE, ECL_INT_TYPE, ECL_DOUBLE_TYPE,
                                            ECL_BOOL_TYPE, ECL_MESS_TYPE};

    static public int ECL_STRING_LENGTH = 8;
    static public int ECL_TYPE_LENGTH = 4;

    static public int BLOCKSIZE_NUMERIC = 1000;
    static public int BLOCKSIZE_CHAR = 105;

    static public boolean isValidKey(String key)
    {
        if(key == null) return false;
        if((key.length() != ECL_STRING_LENGTH))
            return false;
        /*
        for(int i=0; i<ECL_KEYWORDS.length; i++)
        {
            if(key.contains(ECL_KEYWORDS[i]))
                return true;
        }
        System.out.println("Keyword not found....");
        */
        return true;
    }

    class EclipseRecord
    {
        int size = 0;
        int sizeof_ctype = 0;
        byte ecl_type = -1;
        public String header = null;
        byte[] data = null;
        boolean shared = false;
        long offset = -1;

        public EclipseRecord()
        {
        }

        public EclipseRecord(int _size, int _sizeof_ctype, byte _type, byte[] _hdr, byte[] _data, boolean _shared, long _offset)
        {
            size = _size;
            sizeof_ctype = _sizeof_ctype;
            ecl_type = _type;
            header = new String(_hdr);
            data = _data;
            shared = _shared;
            offset = _offset;
        }

        public byte getEclFmt() { return ecl_type; }
        public int getSize() { return size; }
        public long getOffset() { return offset; }
        public boolean getShared() { return shared; }
        public byte[] getData() { return data; }
        public String getHeader() { return header;}
        public int getTypeSize() { return sizeof_ctype; }

        public void setEclFmt(byte _type) { ecl_type = _type; setSizeOfType();}
        public void setSize(int _size) { size = _size; }
        public void setOffset(long _offset) { offset = _offset; }
        public void setShared(boolean _shared) { shared = _shared; }
        public void setData(byte[] _data) { data = _data; }
        public void setHeader(String header) { this.header = header;}
        public void setTypeSize(int size) { sizeof_ctype = size; }

        public byte getEclFmtFromString(String fmt)
        {
            for(int i=0; i<ECL_TYPE_NAMES.length; i++)
            {
                if(ECL_TYPE_NAMES[i].equalsIgnoreCase(fmt))
                    return ECL_TYPES[i];
            }
            return -1;
        }

        public int getBlockSize()
        {
            if (ecl_type == ECL_CHAR_TYPE)
                return BLOCKSIZE_CHAR;
            else if (ecl_type == ECL_MESS_TYPE)
                return BLOCKSIZE_CHAR;
            else
                return BLOCKSIZE_NUMERIC;
        }

        public int setSizeOfType()
        {
            switch (ecl_type)
            {
                case(ECL_CHAR_TYPE):
                    sizeof_ctype = (ECL_STRING_LENGTH + 1); /* One element of character data is a string section of 8 characters + \0. */
                    break;
                case(ECL_FLOAT_TYPE):
                    sizeof_ctype = 4;
                    break;
                case(ECL_DOUBLE_TYPE):
                    sizeof_ctype = 8;
                    break;
                case(ECL_INT_TYPE):
                    sizeof_ctype = 4;
                    break;
                case(ECL_BOOL_TYPE):
                    sizeof_ctype = 4; // The ECL_BOOL_TYPE type is internally implemented as an integer - and not a bool.
                    break;
                case(ECL_MESS_TYPE):
                    sizeof_ctype = 1;
                    break;
                default:
            }
            return sizeof_ctype;
        }

        public String getDataAsString() { return new String(data); }

        public int getIntegerAt(int idx)
        {
            if(idx*getTypeSize() > data.length)
                return -1;
            else
                return StsMath.convertIntBytes(data, idx*getTypeSize());
        }

        public float getFloatAt(int idx)
        {
            if(idx*getTypeSize() > data.length)
                return -1;
            else
                return StsMath.convertIBMFloatBytes(data, idx*getTypeSize());
        }

        public StsFloatVector getFloatVector()
        {
            StsFloatVector vector = new StsFloatVector(getNumValues());
            for(int i=0; i<getNumValues(); i++)
            {
                vector.append(getFloatAt(i));
            }
            return vector;
        }

        public int getNumValues()
        {
            return data.length/getTypeSize();
        }

        public String toString()
        {
            return header;
        }
    }
/*
    public long getTimeClosestTo(long time)
    {
        getTimes();
        if(times == null)
            return -1;
        for(int i=1; i<times.length; i++)
        {
            if(time > times[i-1] && time < times[i])
                return times[i-1];
        }
        return -1;
    }
*/
}