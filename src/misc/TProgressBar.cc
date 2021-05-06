//
// Project     : HLib
// File        : TProgressBar.cc
// Description : class for showing the progress in various formats
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <iostream>
#include <string>

#include <boost/format.hpp>

#include "hpro/base/System.hh"
#include "hpro/base/error.hh"

#include "hpro/misc/TProgressBar.hh"

namespace HLIB
{

using std::string;

////////////////////////////////////////////////
//
// local variables
//
////////////////////////////////////////////////

namespace
{

// default length of bar
const uint DEFAULT_BAR_LEN = 40;

}// namespace anonymous

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//
// TProgressBar
//
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//
// constructor and destructor
//

TProgressBar::TProgressBar ()
        : _min(0), _max(100), _curr(0), 
          _is_init( false ), _do_cancel(false)
{
}

TProgressBar::~TProgressBar ()
{
    finish();
}

//
// initialise status
//
void
TProgressBar::init ( const double amin, const double amax, const double acurr )
{
    TScopedLock  lock( mutex() );
    
    _min  = amin;
    _max  = amax;
    _curr = acurr;

    // consistency checks: _min ≤ _curr ≤ _max
    _max  = std::max( _min, _max );
    _curr = std::min( std::max( _min, _curr ), _max );
    
    _is_init = true;

    update();
}

//
// reset status, e.g. set new values without initialisation
//
void
TProgressBar::reset ( const double amin, const double amax, const double acurr )
{
    TScopedLock  lock( mutex() );
    
    _min  = amin;
    _max  = amax;
    _curr = acurr;
    
    // consistency checks: _min ≤ _curr ≤ _max
    _max  = std::max( _min, _max );
    _curr = std::min( std::max( _min, _curr ), _max );
    
    _is_init = true;
    
    update();
}

//
// advance progress
//
void
TProgressBar::advance ( const double f )
{
    TScopedLock  lock( mutex() );

    if ( f != 0.0 )
    {
        _curr += f;

        // consistency check: _min ≤ _curr ≤ _max
        // _curr = std::min( std::max( _min, _curr ), _max );
        
        update();
    }// if
}

//
// finish output
//
void
TProgressBar::finish ()
{
    TScopedLock  lock( mutex() );
    
    _is_init = false;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//
// TStdProgressBar
//
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

#if defined(LINUX) || defined(SUNOS) || defined(AIX) || \
    defined(TRU64) || defined(HPUX)  || defined(DARWIN)
const string STD_FORMAT = "%b %p %e (%m)";
#else
const string STD_FORMAT = "%b %p %e";
#endif


//
// constructor and destructor
//
TConsoleProgressBar::TConsoleProgressBar ( const term_charset_t  /* acharset */,
                                           const term_color_t    /* acolor */ )
        : _out( & ( std::cout ) )
        , _text_len(0)
        , _old_val(-1)
        , _charset( ascii_charset )
        , _color_mode( no_color )
{
    if ( _out == nullptr )
        HERROR( ERR_ARG, "(TConsoleProgressBar) ctor", "nullptr output stream given" );
    
    _format = STD_FORMAT;

    _indicator_status = 0;
}

TConsoleProgressBar::TConsoleProgressBar ( const string &        format,
                                           const term_charset_t  /* acharset */,
                                           const term_color_t    /* acolor */ )
        : _out( & ( std::cout ) )
        , _text_len(0)
        , _old_val(-1)
        , _charset( ascii_charset )
        , _color_mode( no_color )
{
    if ( _out == nullptr )
        HERROR( ERR_ARG, "(TConsoleProgressBar) ctor", "nullptr output stream given" );
    
    if ( format == "" ) _format = STD_FORMAT;
    else                _format = format;

    _indicator_status = 0;
}

TConsoleProgressBar::TConsoleProgressBar ( std::ostream &        out_stream,
                                           const term_charset_t  /* acharset */,
                                           const term_color_t    /* acolor */ )
        : _out( & out_stream )
        , _text_len(0)
        , _old_val(-1)
        , _charset( ascii_charset )
        , _color_mode( no_color )
{
    if ( _out == nullptr )
        HERROR( ERR_ARG, "(TConsoleProgressBar) ctor", "nullptr output stream given" );
    
    _format = STD_FORMAT;

    _indicator_status = 0;
}

TConsoleProgressBar::TConsoleProgressBar ( std::ostream &        out_stream,
                                           const string &        format,
                                           const term_charset_t  /* acharset */,
                                           const term_color_t    /* acolor */ )
        : _out( & out_stream )
        , _text_len(0)
        , _old_val(-1)
        , _charset( ascii_charset )
        , _color_mode( no_color )
{
    if ( _out == nullptr )
        HERROR( ERR_ARG, "(TConsoleProgressBar) ctor", "nullptr output stream given" );
    
    if ( format == "" ) _format = STD_FORMAT;
    else                _format = format;

    _indicator_status = 0;
}

TConsoleProgressBar::~TConsoleProgressBar ()
{
    finish();
}

//
// initialise status
//
void
TConsoleProgressBar::init ( const double  amin,
                            const double  amax,
                            const double  acurr )
{
    {
        TScopedLock  lock( mutex() );

        _old_val  = -1;
        _text_len = 0;
    }

    TProgressBar::init( amin, amax, acurr );
}

//
// finish output
//
void
TConsoleProgressBar::finish ()
{
    {
        TScopedLock  lock( mutex() );

        clear_line();
    }
    
    TProgressBar::finish();
}

//
// print progress
//
void
TConsoleProgressBar::update ()
{
    const double  percent = 100.0 * percentage();

    // check if update is really necessary
    if ( uint(percent) == uint(_old_val) )
        return;
    
    _old_val = percent;

    //
    // print actual bar by parsing format string and inserting
    // wished values
    //

    const size_t  format_len = _format.length();
    string        output, number;
    size_t        len = 0;
    uint          field_len;
    
    for ( uint i = 0; i < format_len; i++ )
    {
        if ( _format[i] == '%' )
        {
            // check for end of string
            if ( ++i >= format_len )
                break;

            // look for length field
            number = "";
            while (( i < format_len ) && ( isdigit( _format[i] ) ))
                number += _format[i++];

            if ( number != "" ) field_len = atoi( number.c_str() );
            else                field_len = 0;
                
            // again check for end of string
            if ( i >= format_len )
                break;

            if ( _format[i] == 'b' )
            {
                if ( field_len == 0 )
                    field_len = DEFAULT_BAR_LEN;
                
                //
                // insert bar
                //

                const double  bar_len      = double( field_len ) * percent / 100.0;
                const uint    full_bar_len = uint( std::floor( bar_len ) );
                char          empty_char   = '_';

                if ( _color_mode == use_color )
                {
                    empty_char = ' ';
                }// if
                
                if ( _charset == unicode_charset )
                {
                    // const uint  split_ratio  = uint( 3.0 * ( bar_len - double( full_bar_len ) ) );
                    const uint  split_ratio  = uint( 4.0 * ( bar_len - double( full_bar_len ) ) );
                    
                    for ( uint j = 0; j < field_len; j++ )
                    {
                        if ( j < full_bar_len  )
                        {
                            output += "█";
                        }// if
                        else if ( j == full_bar_len )
                        {
                            switch ( split_ratio )
                            {
                                case 0  : output += empty_char; break;
                                case 1  : output += "▎"; break;
                                case 2  : output += "▌"; break;
                                case 3  : output += "▊"; break;
                                default : output += " "; break;

                                // case 0  : output += "░"; break;
                                // case 1  : output += "▒"; break;
                                // case 2  : output += "▓"; break;
                                // default : output += "░"; break;
                            }// switch
                        }// if
                        else
                        {
                            output += empty_char;
                            // output += "░";
                        }// else
                    }// for
                
                    len += field_len;
                }// if
                else if ( _charset == ascii_charset )
                {
                    const uint  split_ratio  = uint( 2.0 * ( bar_len - double( full_bar_len ) ) );
                    
                    output += '[';
                    len++;
                
                    for ( uint j = 0; j < field_len; j++ )
                    {
                        if ( j < full_bar_len  )
                        {
                            output += '=';
                        }// if
                        else if ( j == full_bar_len )
                        {
                            switch ( split_ratio )
                            {
                                case 0  : output += ' '; break;
                                case 1  : output += '-'; break;
                                default : output += ' '; break;
                            }// switch
                        }// if
                        else
                        {
                            output += ' ';
                        }// else
                    }// for
                
                    len += field_len;
                    
                    output += ']';
                    len++;
                }// else
            }// if
            else if ( _format[i] == 'B' )
            {
                //
                // insert single character bar
                //

                if ( _charset == unicode_charset )
                {
                    // height in (stretched) 1/8 steps
                    const uint  height = Math::min( uint(16), uint(double(16) * percent / 100.0) );

                    switch ( height )
                    {
                        case  0 : output += " "; break;
                        case  1 :
                        case  2 : output += "▁"; break;
                        case  3 :
                        case  4 : output += "▂"; break;
                        case  5 :
                        case  6 : output += "▃"; break;
                        case  7 :
                        case  8 : output += "▄"; break;
                        case  9 :
                        case 10 : output += "▅"; break;
                        case 11 :
                        case 12 : output += "▆"; break;
                        case 13 :
                        case 14 : output += "▇"; break;
                        case 15 :
                        case 16 : output += "█"; break;
                    }// switch
                
                    len++;
                }// if
                else if ( _charset == ascii_charset )
                {
                    // height in (stretched) 1/4 steps
                    const uint  height = Math::min( uint(8), uint(double(8) * percent / 100.0) );
                    
                    switch ( height )
                    {
                        case  0 : output += " "; break;
                        case  1 : 
                        case  2 : output += "."; break;
                        case  3 : 
                        case  4 : output += "o"; break;
                        case  5 : 
                        case  6 : output += "O"; break;
                        case  7 : 
                        case  8 : output += "0"; break;
                    }// switch
                
                    len++;
                }// else
            }// if
            else if ( _format[i] == 'i' )
            {
                //
                // show progress indicator
                //

                if ( _charset == unicode_charset )
                {
                    switch ( _indicator_status )
                    {
                        case 0 : output += "◐"; _indicator_status++; break;
                        case 1 : output += "◓"; _indicator_status++; break;
                        case 2 : output += "◑"; _indicator_status++; break;
                        case 3 : output += "◒"; _indicator_status = 0; break;
                    }// switch
                }// if
                else if ( _charset == ascii_charset )
                {
                    switch ( _indicator_status )
                    {
                        case 0 : output += "/";  _indicator_status++; break;
                        case 1 : output += "-";  _indicator_status++; break;
                        case 2 : output += "\\"; _indicator_status++; break;
                        case 3 : output += "|";  _indicator_status = 0; break;
                    }// switch
                }// else

                len++;
            }// if
            else if ( _format[i] == 'p' )
            {
                //
                // insert percentage info
                //
                
                const std::string  s = boost::str( boost::format( "%3d%%" ) % int( percent ) );
                
                output += s;
                len    += s.length();
            }// if
            else if ( _format[i] == 'm' )
            {
                //
                // insert memory
                //
                
                const std::string  s = Mem::to_string();

                output += s;
                len    += s.length();
            }// if
            else if ( _format[i] == '%' )
            {
                // percent sign
                output += '%';
                ++len;
            }// if
        }// if
        else
        {
            output += _format[i];
            ++len;
        }// else
    }// for

    //
    // clear line and write new output
    //
    
    clear_line();
    *_out << output << std::flush;

    _text_len = uint(len);
}

//
// clear current output line
//
void
TConsoleProgressBar::clear_line ()
{
#if 0
    
    for ( uint  i = 0; i < _text_len; ++i )
        *_out << "\b \b";
    
#else

    *_out << '\r';
    
    for ( uint  i = 0; i < _text_len; ++i )
        *_out << ' ';
    
    *_out << '\r' << std::flush;

#endif
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//
// TPartProgressBar
//
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//
// ctor
//
TPartProgressBar::TPartProgressBar ( TProgressBar  *  parent,
                                     const double     tot_part )
        : _parent( parent )
        , _total_part( tot_part )
{
    if ( _parent == nullptr )
        HERROR( ERR_ARG, "(TPartProgressBar) ctor", "parent is NULL" );
}

//
// advance progress by \a f
//
void
TPartProgressBar::advance ( const double  f )
{
    TScopedLock  lock( mutex() );
    
    if ( f != 0.0 )
    {
        _curr += f;
        
        _parent->advance( _total_part * ( f / ( _max - _min ) ) );
    }// if
}

//
// react upon change in status
//
void
TPartProgressBar::update ()
{
    // no update
}

}// namespace HLIB
