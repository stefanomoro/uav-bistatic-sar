#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: Not titled yet
# Author: nuc1
# GNU Radio version: v3.8.5.0-5-g982205bd

from gnuradio import blocks
from gnuradio import gr
from gnuradio.filter import firdes
import sys
import signal
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio import eng_notation
from gnuradio import uhd
import time
from datetime import datetime


class usrp_rx_working(gr.top_block):

    def __init__(self, freq=2e9, samp_rate=56e6, file_sink_folder='/mnt/ramdisk/'):
        gr.top_block.__init__(self, "USRP RX for SAR")

        ##################################################
        # Parameters
        ##################################################
        self.freq = freq
        self.samp_rate = samp_rate

        ##################################################
        # Blocks
        ##################################################
        self.uhd_usrp_source_0 = uhd.usrp_source(
            ",".join(("", 'num_recv_frames=1024')),
            uhd.stream_args(
                cpu_format="sc16",
                args='',
                channels=list(range(0, 1)),
            ),
        )
        self.uhd_usrp_source_0.set_time_source('gpsdo', 0)
        self.uhd_usrp_source_0.set_clock_source('gpsdo', 0)
        self.uhd_usrp_source_0.set_center_freq(freq, 0)
        self.uhd_usrp_source_0.set_rx_agc(False, 0)
        self.uhd_usrp_source_0.set_gain(60, 0)
        self.uhd_usrp_source_0.set_antenna('TX/RX', 0)
        self.uhd_usrp_source_0.set_bandwidth(samp_rate, 0)
        self.uhd_usrp_source_0.set_samp_rate(samp_rate)
        self.uhd_usrp_source_0.set_time_unknown_pps(uhd.time_spec())
        self.uhd_usrp_source_0.set_min_output_buffer(100000000)
        self.uhd_usrp_source_0.set_max_output_buffer(8192)
        
        gps_time = self.uhd_usrp_source_0.get_mboard_sensor("gps_time", 0)
        #time_at_last_pps = self.uhd_usrp_source_0.get_time_now().get_tick_count()#get_real_secs()
        print("TIME {}".format(gps_time))
        
        #if not(gps_time.to_int()):
        file_name = file_sink_folder + getFileName(self.freq,self.samp_rate)
        #else:
            #file_name = file_sink_folder + str(gps_time.to_int()) + ".dat"
        
        
        self.blocks_stream_to_vector_0 = blocks.stream_to_vector(
            gr.sizeof_short*2, 1)
        self.blocks_file_sink_1 = blocks.file_sink(
            gr.sizeof_short*2, file_name, False)
        self.blocks_file_sink_1.set_unbuffered(False)

        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_stream_to_vector_0, 0),
                     (self.blocks_file_sink_1, 0))
        self.connect((self.uhd_usrp_source_0, 0),
                     (self.blocks_stream_to_vector_0, 0))

    def get_freq(self):
        return self.freq

    def set_freq(self, freq):
        self.freq = freq
        self.uhd_usrp_source_0.set_center_freq(self.freq, 0)

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.uhd_usrp_source_0.set_samp_rate(self.samp_rate)
        self.uhd_usrp_source_0.set_bandwidth(self.samp_rate, 0)


def argument_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "--freq", dest="freq", type=eng_float, default="2.0G",
        help="Set freq [default=%(default)r]")
    parser.add_argument(
        "--samp-rate", dest="samp_rate", type=eng_float, default="56.0M",
        help="Set samp_rate [default=%(default)r]")
    return parser


def main(top_block_cls=usrp_rx_working, options=None):
    if options is None:
        options = argument_parser().parse_args()
    #file_name = getFileName(freq=options.freq, samp_rate=options.samp_rate)
    tb = top_block_cls(freq=options.freq,
                       samp_rate=options.samp_rate)

    def sig_handler(sig=None, frame=None):
        tb.stop()
        tb.wait()

        sys.exit(0)

    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    tb.start()

    try:
        input('Press Enter to quit: ')
    except EOFError:
        pass
    tb.stop()
    tb.wait()


def getFileName(freq, samp_rate):
    fname = datetime.now().strftime("%Y-%m-%d %H:%M:%S") + \
            "_f" + str(int(freq/1e6)) + "_s"+str(int(samp_rate/1e6)) + ".dat"

    return fname


if __name__ == '__main__':
    main()
