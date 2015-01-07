package com.Sts.IO;

import java.io.*;



/**
	Buffers access to {@link java.io.RandomAccessFile}, which enormously speeds I/O.
	Also consider Java 1.4's {@link java.nio}.
	This is class is useful for old code and with libraries,
	as one can simply make a <code>new BuffeadedRandomAccessFile()</code> instance and pass to existing code wherever a RandomAccessFile is expected,
	without rewriting that code.
	This class has been used extensively and tested extensively with long sequence of random file operations.
	Nevertheless: NO WARRENTY.

	<p><em>WARNING!</em>  This class is not a complete drop-in replacement for RandomAccessFile.
	Do not use the two methods inherited from {@link RandomAccessFile}:
		{@link RandomAccessFile#writeBytes(String s)} and
		{@link RandomAccessFile#writeChars(String s)}.
		They are not compatible with buffered use.
	<em>Instead use</em> the methods {@link #writeString8(String)} and {@link #writeString16(String)}, respectively.
	(Also, {@link java.lang.String}s can also be written as UTF-8 with {@link RandomAccessFile#writeUTF(String)}).

	<p>The technical reason for the warning is that the methods themselves are <code>final</code>
	and ultimately write with the <code>private</code>, and hence non-overrideable, method <code>writeBytes(byte[] b, int off, int len)</code>.
	Other methods all other write methods use <code>write()</code>/<code>write(byte[])/write(byte[],int,int)</code> for low-level byte writing,
	and since those methods are overriden in this class, all is well.
	This class could be written with the exact same interface as {@link java.io.RandomAccessFile}
	(by nesting rather than subclassing).
	However this would require all code using the class to be recompiled, which is often not an option, especially in third-party code.
	A small change by Sun would make those two methods valid in this class as well --
	<a href='http://developer.java.sun.com/developer/bugParade/bugs/4400941.html'>vote for it</a>.

<!--
	maybe add peek()
	setting for back buffer, so can read file backwards efficiently?
	when seek pos < getFilePointer, maybe seek 100 bytes or so before so can efficiently seek backwards
	X constructor that wraps existing RandomAccessFile => would means switching from inheriting functionality to composition
-->

	@author T.A. Phelps
	@version $Revision: 1.5 $ $Date: 2010/09/22 17:26:12 $
*/
public class StsNioRandomAccessFile extends RandomAccessFile
{
    static final int DEFAULT_BUFSIZE = 8*1024;    // adjust to most common disk block size
    static final int MIN_READ = 8*1024;

    /*private or protected*/byte[] buf_;
    /** Current point within buffer. */
    private int bufi_;
    /** File pointer of start of buffer, so fp_ + bufi_ == {@link #getFilePointer()}. */
    private long fp_ = 0;
    /** Length of valid portion of buffer, buf_[0..limit_-1]. */
    private int limit_ = 0;
    /** Does buffer have changed but unwritten data? */
    private boolean dirty_ = false;
    /** {@link RandomAccessFile#length()} also slow, so cache it.  Updated on writes. */
    private long length_=-1;



    public StsNioRandomAccessFile(File file, String mode) throws FileNotFoundException { this(file, mode, DEFAULT_BUFSIZE); }
    public StsNioRandomAccessFile(String name, String mode) throws FileNotFoundException { this(name, mode, DEFAULT_BUFSIZE); }
    public StsNioRandomAccessFile(String name, String mode, int bufsize) throws FileNotFoundException { this(new File(name), mode, bufsize); }

    public StsNioRandomAccessFile(File file, String mode, int bufsize) throws FileNotFoundException
    {
        super(file, mode);

        int maxbuf = Math.max(bufsize,1);
        if ("r".equals(mode)) try { maxbuf = (int)Math.min(maxbuf, length()); } catch (IOException failok) {}
        buf_ = new byte[maxbuf];
        //System.out.println("bufsize = "+buf_.length);
        try { invalidate(); } catch (IOException canthappen) {}
    }


    /** Flush and update fp_, then invalidate buffer contents to force buffer refresh on next read/write. */
    /*protected? public?*/ void invalidate() throws IOException
    {
        flush();
        fp_ += bufi_;   // set fp_ to file pointer, because bufi_ going to 0, so fp_ + bufi_ wrong
        bufi_ = limit_ = 0;
        //dirty_=false; // done in flush()
    }

    public int read(byte b[]) throws IOException { return read(b, 0, b.length); }     // not the default implementation that you'd think, but shouldn't realy on this anyhow

    public int read(byte b[], int off, int len) throws IOException
    {
//	    assert b!=null && off>=0 && len>=0 && len+off <= b.length;
//	    assert limit_>=bufi_: limit_+" < "+bufi_;

        int got;
        int already = Math.min(limit_ - bufi_, len);
        if (already > 0)
        {
            System.arraycopy(buf_,bufi_, b,off, already);
            bufi_ += already;
            got = already;

        }
        else if (len == 0) got = 0;
        else if (limit_ == -1) got = -1;

        else
        {
//            assert bufi_ == limit_;
            invalidate();   // could switch from write to read
            limit_ = super.read(buf_, 0, buf_.length);
            got = read(b, off, len);    // recurse
        }

//        assert got==-1 || got>=1 || len==0;   // "blocks until at least one byte of input is available"
        return got;
    }

    public int read() throws IOException {
        int b;

        if (bufi_ < limit_) b = (buf_[bufi_++] & 0xff);
        else if (limit_ == -1) b = -1;
        else
        {
//            assert bufi_ == limit_;
            invalidate();
            limit_ = super.read(buf_, 0, buf_.length);
            b = read(); // recurse
        }
        return b;
    }




    // X private void writeBytes(byte[] b, int off, int len) throws IOException => would catch all writes, but can't override a private!

    /**
    Write a {@link java.lang.String} as a sequence of one byte per character, discarding the high byte.
    Like {@link #writeBytes(String)}, which for technical reasons must not be used.
    */
    public void writeString8(String s) throws IOException
    {
        // DataOutputStream.write(String) does it this way:
        /*if (s!=null)*/
        for (int i=0,imax=s.length(); i<imax; i++) write(s.charAt(i));     // no alloc, no deprecated method
    }

    /**
    Same as {@link RandomAccessFile#writeChar(int)} for each character of the String in sequence.
    Like {@link RandomAccessFile#writeChars(String)}, which for technical reasons must not be used.
    */
    public void writeString16(String s) throws IOException {
        /*if (s!=null)*/
        for (int i=0,imax=s.length(); i<imax; i++) writeChar(s.charAt(i));
        // RandomAccessFile.writeChars(String) builds byte array and writes that as block, but not so important with buffered i/o
    }

    public void write(byte b[]) throws IOException { write(b, 0, b.length); }     // identical to inherited but be safe

    public void write(byte[] b, int off, int len) throws IOException
    {
//        assert b!=null && off>=0 && len>=0 && off+len <= b.length;
        if (len == 0) return;   // degenerate case

        // write what you can to buffer
        int space = Math.min(buf_.length - bufi_, len);
        if (space > 0)
        {
            System.arraycopy(b,off, buf_,bufi_, space);
            bufi_ += space;
            limit_ = Math.max(bufi_, limit_);   // not necessarily at end of buffer
            dirty_ = true;

            off += space; len -= space;
        }

        if (len > 0)
        {
//            assert bufi_ == limit_;
            invalidate();   // if writing into middle of dirty buf, can uselessly write bufi_ .. limit_, but that's not much and only if mix write() and write(byte[])

            if (len < buf_.length) write(b,off, len);   // recurse
            else
            {
//                assert fp_ == super.getFilePointer(): fp_+" != "+super.getFilePointer();
                //super.seek(fp_);    // can seek back in file -- but synced up with invalidate()
                super.write(b, off, len);  // can't fit into buffer no how
                fp_ += len;
            }
        }
        if (fp_ + bufi_ > length_) length_ = fp_ + bufi_;
    }

    public void write(int b) throws IOException
    {
        if (bufi_ < limit_) {}   // somewhere in already written space
        else if (limit_ < buf_.length)
        {
            // fits in unused part of buffer
            if (limit_==-1) limit_=0;   // working out of EOF hole from read()
//            assert bufi_==limit_;
            limit_++;
        }
        else
        {    // reset buffer
//            assert bufi_==limit_;
            invalidate();
            limit_=1;   // 1 after write byte below
        }

        buf_[bufi_++] = (byte)b;

        dirty_ = true;
        if (fp_ + bufi_ > length_) length_ = fp_ + bufi_;
    }

    //public void writeBytes(String s) {} -- final!

    /** Flush unwritten data from buffer to file.  Buffer and file pointer are unchanged. */
    public void flush() throws IOException
    {
        if (!dirty_) return;

        // could (1) write only at end of file or (2) read, seek, write, then read more which forces flush
        /*if (bufi_ > 0)--almost always*/
        super.seek(fp_);     // generally super at fp_+limit_

        super.write(buf_, 0, limit_);   // usually only 0..bufi_ valid, but could have seek() within buf

        //fp_ += bufi_; => doesn't affect buffer
        //invalidate(); => buffer still valid
        dirty_=false;
    }

    public long getFilePointer() throws IOException
    {
        return fp_ + (limit_>=0? bufi_: 0);     // limit_==-1 at EOF
    }

    // LATER: guarantee skips n bytes, except if can't because at end of file -- but I suppose read(byte[]) craps out in same cases as skipBytes()
    public int skipBytes(int n) throws IOException
    {
        if (n<=0) return 0;  // in RandomAccessFile contract
        n = (int)Math.min(length() - getFilePointer(), n);

        int cnt = Math.min(limit_ - bufi_, n);  // try to satisfy in current buffer
        bufi_ += cnt; n -= cnt;
        if (n > 0)
        {
            invalidate();
            cnt += super.skipBytes(n);
        //fp_ += ss; => does a seek internally so already updated!
        }
        return cnt;
    }

    public void setLength(long newLength) throws IOException
    {
//        assert newLength >=0: newLength;
        if (newLength == length()) return;  // no change

        invalidate();   // all bets are off
        super.setLength(newLength);
        length_ = newLength;
        if (fp_ > newLength) fp_ = newLength;   // file position truncated to EOF on Windows 98, though Java API doesn't define
        super.seek(fp_);
        //fp_ = super.getFilePointer(); -- wrong because out of sync
    }

    public long length() throws IOException
    {
        if (length_ < 0) length_ = super.length();  // HUGELY slow on Windows 98
        return length_;
    }

    public void seek(long pos) throws IOException
    {
        long end = length();  //assert pos <= length(): pos+" "+length(); => can seek beyond end

        if (fp_ <= pos && pos < fp_+limit_)
        {   // within current buffer
            bufi_ = (int)(pos - fp_);

        }
        else if (end - MIN_READ <= pos && pos <= end)
        {
            // special case: buffer at leasat MIN_READ, so can efficiently seek backwards from EOF as for DVI and PDF.  full buffer size can be wasteful if havae 2MB buf and just read a few bytes at end of file
            //always cache up to buffer size,
            invalidate();
            int max = Math.min((int)Math.min(buf_.length, end), MIN_READ);  // readFully() zaps limit_ setting, so use local here
            super.seek(end - max);
            //System.arraycopy(); // LATER: if overlap, reuse part of buffer
            readFully(buf_, 0, max);     // fill here to maintain file pointer -- could do this at buffer fill time, but this the more unusual case (seeking near EOF)
            limit_ = max; fp_ = end - limit_; bufi_ = (int)(pos - fp_);

        }
        else
        {    // fault in, or seek past end of file
            invalidate();
            super.seek(pos);    // previous 100 bytes too so can scan backwards efficiently
            fp_ = pos;
        }
//        assert fp_ + bufi_ == pos && bufi_>=0: fp_+" "+bufi_+" "+pos;
    }

    public void close() throws IOException
    {
        flush();
        super.close();
    }
}