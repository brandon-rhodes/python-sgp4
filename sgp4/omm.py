"""Support for the new Orbit Mean-Elements Message format for TLE data."""

import csv
import xml.etree.ElementTree as ET
from datetime import datetime
from math import pi
from sgp4.api import WGS72

def parse_csv(file):
    return csv.DictReader(file)

def parse_xml(file):
    root = ET.parse(file).getroot()
    for segment in root.findall('.//segment'):
        metadata = segment.find('metadata')
        data = segment.find('data')
        meanElements = data.find('meanElements')
        tleParameters = data.find('tleParameters')
        fields = {}
        for element in metadata, meanElements, tleParameters:
            fields.update((field.tag, field.text) for field in element)
        yield fields

_epoch0 = datetime(1949, 12, 31)
_to_radians = pi / 180.0
_ndot_units = 1036800.0 / pi  # See SGP4.cpp for details.
_nddot_units = 2985984000.0 / 2.0 / pi  # See SGP4.cpp for details.

def initialize(sat, fields):
    sat.classification = fields['CLASSIFICATION_TYPE']
    sat.intldesg = fields['OBJECT_ID'][2:].replace('-', '')
    sat.ephtype = int(fields['EPHEMERIS_TYPE'])
    sat.elnum = int(fields['ELEMENT_SET_NO'])
    sat.revnum = int(fields['REV_AT_EPOCH'])

    epoch_datetime = datetime.strptime(fields['EPOCH'], '%Y-%m-%dT%H:%M:%S.%f')
    epoch = (epoch_datetime - _epoch0).total_seconds() / 86400.0

    argpo = float(fields['ARG_OF_PERICENTER']) * _to_radians
    bstar = float(fields['BSTAR'])
    ecco = float(fields['ECCENTRICITY'])
    inclo = float(fields['INCLINATION']) * _to_radians
    mo = float(fields['MEAN_ANOMALY']) * _to_radians
    nddot = float(fields['MEAN_MOTION_DDOT']) / _nddot_units
    ndot = float(fields['MEAN_MOTION_DOT']) / _ndot_units
    no_kozai = float(fields['MEAN_MOTION']) / 720.0 * pi
    nodeo = float(fields['RA_OF_ASC_NODE']) * _to_radians
    satnum = int(fields['NORAD_CAT_ID'])

    sat.sgp4init(WGS72, 'i', satnum, epoch, bstar, ndot, nddot, ecco,
                 argpo, inclo, mo, no_kozai, nodeo)
