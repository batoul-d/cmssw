/**\class CSCWireDigi
 *
 * Digi for CSC anode wires.
 *
 */

#include "DataFormats/CSCDigi/interface/CSCWireDigi.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iostream>
#include <cstdint>

  /// Constructors

CSCWireDigi::CSCWireDigi (int wire, unsigned int tbinb){
  /// Wire group number in the digi (16 lower bits from the wire group number)
  wire_  = wire & 0x0000FFFF;
  tbinb_ = tbinb;
  /// BX in the wire digis (16 upper bits from the wire group number)
  wireBX_= (wire >> 16) & 0x0000FFFF;
  /// BX wire group combination 
  wireBXandWires_=wire;
}

   /// Default
CSCWireDigi::CSCWireDigi (){
  wire_ = 0;
  tbinb_ = 0;
  wireBX_=0;
  wireBXandWires_=0;
}

  /// return tbin number (obsolete, use getTimeBin() instead)
int CSCWireDigi::getBeamCrossingTag() const {
  return getTimeBin();
}
  /// return first tbin ON number
int CSCWireDigi::getTimeBin() const {
  uint32_t tbit=1;
  int tbin=-1;
  for(int i=0;i<32;i++) {
    if(tbit & tbinb_) tbin=i;
    if(tbin>-1) break;
    tbit=tbit<<1;
  }
  return tbin;
}
  /// return vector of time bins ON
std::vector<int> CSCWireDigi::getTimeBinsOn() const {
  std::vector<int> tbins;
  uint32_t tbit=tbinb_;
  uint32_t one=1;
  for(int i=0;i<32;i++) {
    if(tbit & one) tbins.push_back(i);
    tbit=tbit>>1;
    if(tbit==0) break;
  }
  return tbins;
}

  /// Debug

void CSCWireDigi::print() const {
  std::ostringstream ost;
  ost << "CSCWireDigi | wg " << getWireGroup() << " | "
      << " BX # " << getWireGroupBX() << " | "
      << " BX + Wire " << std::hex << getBXandWireGroup() << " | " << std::dec
      << " First Time Bin On " << getTimeBin() << " | Time Bins On ";
  std::vector<int> tbins=getTimeBinsOn();
  for(unsigned int i=0; i<tbins.size();i++) {ost << tbins[i] << " ";}
  edm::LogVerbatim("CSCDigi") << ost.str();
}

std::ostream & operator<<(std::ostream & o, const CSCWireDigi& digi) {
  o << " CSC Wire " << digi.getWireGroup()
           << " CSC Wire First Time Bin On " << digi.getTimeBin()
           << " CSC Time Bins On ";
  for (unsigned int i = 0; i<digi.getTimeBinsOn().size(); ++i ){
    o <<" " <<digi.getTimeBinsOn()[i]; }
  return o;
}

